"""Export a research-article reader; figures and captions remain one object."""
from pathlib import Path
import re
import nbformat
from nbconvert import HTMLExporter
from bs4 import BeautifulSoup
ROOT=Path(__file__).resolve().parent


def main():
    book=nbformat.read(ROOT/'铜热膨胀研究.ipynb',as_version=4)
    exporter=HTMLExporter(template_name='lab')
    exporter.exclude_input_prompt=True
    exporter.exclude_output_prompt=True
    html,_=exporter.from_notebook_node(book)
    soup=BeautifulSoup(html,'html.parser')
    soup.body['class']=soup.body.get('class',[])+['paper-reader']
    style=soup.new_tag('style')
    style.string=(ROOT/'docs/paper_style.css').read_text(encoding='utf8')
    soup.head.append(style)
    soup.title.string='Thermal Expansion of Single-Crystal Copper | Research Article'
    rendered=soup.select('.jp-Cell')
    assert len(rendered)==len(book.cells),(len(rendered),len(book.cells))
    for index,(cell,element) in enumerate(zip(book.cells,rendered)):
        if cell.metadata.get('paper_figure'):
            number=cell.metadata.paper_figure
            element['id']=f'fig{number}'
            element['class']=element.get('class',[])+['paper-figure']
            content=element.select_one('.jp-RenderedHTMLCommon')
            paragraphs=content.find_all('p',recursive=False)
            image_paragraph=next(p for p in paragraphs if p.find('img'))
            caption_paragraphs=[p for p in paragraphs if p!=image_paragraph]
            figure=soup.new_tag('figure')
            image=image_paragraph.find('img').extract()
            figure.append(image)
            caption=soup.new_tag('figcaption')
            for p in caption_paragraphs:
                caption.append(p.extract())
            figure.append(caption)
            content.clear()
            content.append(figure)
        if cell.cell_type=='code':
            wrapper=element.select_one('.jp-Cell-inputWrapper')
            if wrapper:
                details=soup.new_tag('details',attrs={'class':'paper-code'})
                wrapper.wrap(details)
                summary=soup.new_tag('summary')
                summary.string='Executable analysis · click to expand'
                details.insert(0,summary)
    # A restrained abstract block and bibliography sizing, without faux branding.
    for heading in soup.select('h2'):
        title=heading.get_text().replace('¶','').strip()
        if title=='Abstract':
            heading.find_parent(class_='jp-Cell')['class'].append('paper-abstract')
    for p in soup.select('.jp-RenderedHTMLCommon p'):
        text=p.get_text().strip()
        if re.match(r'^Table\s+\d',text):
            p['class']=p.get('class',[])+['paper-table-caption']
        elif re.match(r'^\[\d+\]',text):
            p['class']=p.get('class',[])+['paper-reference']
            p['id']='ref'+re.match(r'^\[(\d+)\]',text)[1]
        # Link figure references in prose, preserving nested formatting and math.
        if not p.find_parent('figcaption'):
            for node in list(p.find_all(string=re.compile(r'Figure [1-6]\b'))):
                if node.find_parent('a'):
                    continue
                pieces=re.split(r'(Figure [1-6])\b',str(node))
                for piece in pieces:
                    if re.fullmatch(r'Figure [1-6]',piece):
                        link=soup.new_tag('a',href='#fig'+piece[-1])
                        link.string=piece
                        node.insert_before(link)
                    else:
                        node.insert_before(piece)
                node.extract()
    out=ROOT/'results/copper_thermal_expansion_paper.html'
    out.write_text(str(soup),encoding='utf8')
    assert len(soup.select('figure'))==6
    print('Article reader exported with six integrated figures:',out)


if __name__=='__main__':
    main()
