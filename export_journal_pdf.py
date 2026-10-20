"""Typeset the reviewed manuscript as a two-column research reading edition."""
from pathlib import Path
import re, json, subprocess, shutil
import mistune

ROOT=Path(__file__).resolve().parent
PARSER=mistune.create_markdown(renderer='ast',plugins=['table','math'])
SUPER={'⁰':'0','¹':'1','²':'2','³':'3','⁴':'4','⁵':'5','⁶':'6','⁷':'7','⁸':'8','⁹':'9','⁻':'-','⁺':'+'}


def escape(text):
    pieces=re.split(r'([⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺]+)',text)
    out=[]
    for part in pieces:
        if part and all(c in SUPER for c in part):
            out.append(r'\textsuperscript{'+''.join(SUPER[c] for c in part)+'}')
        else:
            out.append(''.join({'\\':r'\textbackslash{}','&':r'\&','%':r'\%','$':r'\$',
                '#':r'\#','_':r'\_','{':r'\{','}':r'\}','~':r'\textasciitilde{}',
                '^':r'\textasciicircum{}'}.get(c,c) for c in part))
    return ''.join(out)


def inline(tokens):
    result=[]
    for token in tokens:
        kind=token['type']
        child=token.get('children',[])
        if kind=='text': result.append(escape(token['raw']))
        elif kind=='strong': result.append(r'\textbf{'+inline(child)+'}')
        elif kind=='emphasis': result.append(r'\emph{'+inline(child)+'}')
        elif kind=='codespan': result.append(r'\texttt{'+escape(token['raw'])+'}')
        elif kind=='inline_math': result.append('$'+token['raw']+'$')
        elif kind=='link':
            url=token['attrs']['url'].replace('%',r'\%').replace('#',r'\#')
            result.append(r'\href{'+url+'}{'+inline(child)+'}')
        elif kind in ['softbreak','linebreak']: result.append(' ')
        elif kind=='inline_html': result.append(' ' if token['raw']=='<br>' else '')
        else: raise ValueError(('Inline token unsupported',token))
    text=''.join(result)
    def citation(match):
        numbers=[i.strip() for i in match[1].split(',')]
        if all(number in {str(n) for n in range(1,13)} for number in numbers):
            return r'\cite{'+','.join('r'+number for number in numbers)+'}'
        return match[0]  # Preserve crystallographic directions such as [100].
    return re.sub(r'\[(\d+(?:,\s*\d+)*)\]',citation,text)


def paragraph(markdown):
    tokens=PARSER(markdown)
    return ' '.join(inline(t.get('children',[])) for t in tokens if t['type']=='paragraph')


def strip_caption(markdown,kind,number):
    return re.sub(r'^\*\*'+kind+r' '+str(number)+r'\.\s*(.*?)\*\*',r'\textbf{\1}',markdown,count=1)


def main():
    source=(ROOT/'docs/manuscript.md').read_text(encoding='utf8').split('## Appendix A.')[0]
    source,references=source.split('## References',1)
    title=source.splitlines()[0].removeprefix('# ')
    abstract_part=source.split('## Abstract',1)[1].split('## 1. Introduction',1)[0]
    abstract,keywords=abstract_part.split('**Keywords:**',1)
    body='## 1. Introduction'+source.split('## 1. Introduction',1)[1]
    body=re.sub(r'\$\$(.*?)\$\$',lambda m:'\n$$\n'+m[1].strip()+'\n$$\n',body,flags=re.S)
    figures={str(f['number']):f for f in json.loads((ROOT/'docs/figure_manifest.json').read_text(encoding='utf8'))}
    blocks=[]
    pending_table=None
    for token in PARSER(body):
        kind=token['type']
        if kind=='blank_line':continue
        if kind=='heading':
            heading=inline(token['children'])
            command='section' if token['attrs']['level']==2 else 'subsection'
            blocks.append('\\'+command+'*{'+heading+'}\n')
        elif kind=='paragraph':
            value=inline(token['children'])
            match=re.match(r'\\textbf\{Table (\d+)\. (.*?)\}(.*)',value,re.S)
            if match:
                pending_table=(match[1],r'\textbf{'+match[2]+'}'+match[3])
            else:
                blocks.append(value+'\n\n')
        elif kind=='block_math':
            blocks.append('\\begin{equation}\n'+token['raw']+'\n\\end{equation}\n')
        elif kind=='block_html':
            match=re.search(r'<!-- figure:(\d+) -->',token['raw'])
            if not match:raise ValueError(token)
            f=figures[match[1]]
            caption=re.sub(r'^\*\*Figure \d+\.\s*(.*?)\*\*',r'**\1**',f['caption'])
            path='../'+(f['pdf'] or f['image'])
            blocks.append('\\begin{figure*}[!tp]\n\\centering\n'
                '\\includegraphics[width=\\textwidth,height=.64\\textheight,keepaspectratio]{'+path+'}\n'
                '\\caption{'+paragraph(caption)+'}\\label{fig'+str(f['number'])+'}\n\\end{figure*}\n')
        elif kind=='table':
            assert pending_table
            number,caption=pending_table
            head=token['children'][0]['children']
            rows=token['children'][1]['children']
            widths=[.23,.38,.39] if len(head)==3 else [.32,.22,.22,.24]
            columns=''.join(r'>{\raggedright\arraybackslash}p{\dimexpr '+f'{w:.3f}'+r'\textwidth-2\tabcolsep\relax}' for w in widths)
            table='\\begin{table*}[!tp]\n\\caption{'+caption+'}\\label{tab'+number+'}\n\\centering\n{\\small\n'
            table+='\\begin{tabular}{'+columns+'}\n\\toprule\n'
            table+=' & '.join(r'\textbf{'+inline(c['children'])+'}' for c in head)+r' \\'+'\n\\midrule\n'
            for row in rows:
                table+=' & '.join(inline(c['children']) for c in row['children'])+r' \\'+'\n'
            table+='\\bottomrule\n\\end{tabular}\n}\n\\end{table*}\n'
            blocks.append(table)
            pending_table=None
        else:raise ValueError(('Unsupported block',token))
    bibliography=[]
    for entry in re.split(r'\n\s*\n',references.strip()):
        match=re.match(r'\[(\d+)\]\s*(.*)',entry,re.S)
        if match:
            bibliography.append(r'\bibitem{r'+match[1]+'} '+paragraph(match[2]))
    assert len(bibliography)==12
    preamble=r'''\documentclass[10pt,a4paper,twocolumn]{article}
\usepackage[margin=18mm,columnsep=6mm,headheight=12pt]{geometry}
\usepackage{fontspec}
\setmainfont{Times New Roman}
\setsansfont{Arial}
\setmonofont{Consolas}
\usepackage{amsmath,amssymb,graphicx,booktabs,array,tabularx}
\usepackage{caption,titlesec,fancyhdr,dblfloatfix,placeins,xurl}
\usepackage[unicode,hidelinks,breaklinks]{hyperref}
\hypersetup{pdftitle={Thermal Expansion of Single-Crystal Copper},pdfauthor={Authors to be completed}}
\DeclareCaptionFont{articlecaption}{\fontsize{8.5}{10.2}\selectfont}
\captionsetup{font={articlecaption},labelfont=bf,labelsep=period,justification=justified,singlelinecheck=false}
\captionsetup[table]{position=top,skip=6pt}
\captionsetup[figure]{skip=7pt}
\titlespacing*{\section}{0pt}{10pt}{5pt}
\titlespacing*{\subsection}{0pt}{8pt}{4pt}
\titleformat{\section}{\bfseries\fontsize{11}{13}\selectfont}{}{0pt}{}
\titleformat{\subsection}{\bfseries\fontsize{10}{12}\selectfont}{}{0pt}{}
\setlength{\parindent}{1em}
\setlength{\parskip}{0pt}
\setlength{\textfloatsep}{12pt plus 2pt minus 2pt}
\setlength{\dbltextfloatsep}{12pt plus 2pt minus 2pt}
\renewcommand{\baselinestretch}{1.0}
\AtBeginDocument{\fontsize{9.5}{11.3}\selectfont}
\renewcommand{\topfraction}{.94}
\renewcommand{\dbltopfraction}{.94}
\renewcommand{\textfraction}{.04}
\renewcommand{\floatpagefraction}{.65}
\renewcommand{\dblfloatpagefraction}{.65}
\setcounter{topnumber}{3}
\setcounter{dbltopnumber}{2}
\setcounter{totalnumber}{4}
\setlength{\emergencystretch}{1.5em}
\makeatletter
\setlength{\@fptop}{0pt}
\setlength{\@dblfptop}{0pt}
\setlength{\@fpsep}{12pt}
\setlength{\@dblfpsep}{12pt}
\makeatother
\urlstyle{same}
\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\footnotesize Copper thermal expansion}
\fancyhead[R]{\footnotesize Research article}
\fancyfoot[C]{\small\thepage}
\renewcommand{\headrulewidth}{.3pt}
\begin{document}
\twocolumn[\begin{@twocolumnfalse}
\begin{center}
{\fontsize{19}{22}\selectfont\bfseries '''+escape(title)+r'''\par}
\vspace{9pt}
{\normalsize [Full names of all group members]\par}
\vspace{4pt}
{\small\itshape [University and School of Materials Science and Engineering]\par}
\end{center}
\vspace{4pt}\hrule\vspace{7pt}
\noindent\textbf{Abstract}\par\smallskip
{\small\noindent '''+paragraph(abstract)+r'''\par}
\vspace{5pt}
{\small\noindent\textbf{Keywords:} '''+paragraph(keywords)+r'''\par}
\vspace{7pt}\hrule\vspace{14pt}
\end{@twocolumnfalse}]
'''
    tex=preamble+'\n'.join(blocks)+'\n\\FloatBarrier\n'
    tex+='\\begin{thebibliography}{99}\n\\fontsize{8.5}{10.2}\\selectfont\n'+'\n\n'.join(bibliography)+'\n\\end{thebibliography}\n\\end{document}\n'
    project=ROOT/'paper_layout'
    project.mkdir(exist_ok=True)
    (project/'main.tex').write_text(tex,encoding='utf8')
    build=ROOT/'.tmp/journal_build'
    build.mkdir(parents=True,exist_ok=True)
    command=['D:/TEX/texlive/2024/bin/windows/xelatex.exe','-interaction=nonstopmode','-halt-on-error',
             '-file-line-error','-output-directory='+str(build),'main.tex']
    for attempt in range(2):
        result=subprocess.run(command,cwd=project,capture_output=True,text=True,encoding='utf8',errors='replace')
        (build/f'compile_{attempt+1}.txt').write_text(result.stdout+result.stderr,encoding='utf8')
        if result.returncode:
            print((result.stdout+result.stderr)[-7000:])
            raise RuntimeError('Journal PDF compilation failed; inspect .tmp/journal_build logs.')
    out=ROOT/'results/copper_thermal_expansion_journal.pdf'
    shutil.copy2(build/'main.pdf',out)
    log=(build/'main.log').read_text(encoding='utf8',errors='replace')
    problems=[line for line in log.splitlines() if 'Overfull' in line or 'Missing character' in line or 'undefined' in line]
    (ROOT/'results/journal_typesetting_check.json').write_text(json.dumps({'compiler':'XeLaTeX',
        'layout':'A4; two columns; 18 mm margins; Times New Roman; full-width figure floats',
        'figures':6,'tables':2,'references':12,'warnings':problems},indent=2),encoding='utf8')
    print('Journal PDF exported:',out)
    print('Typesetting warnings:',problems)


if __name__=='__main__':
    main()
