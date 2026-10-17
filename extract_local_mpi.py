"""Extract Microsoft-signed redistributable payload to project-local runtime.

Requires the official msmpisetup.exe already downloaded to .tmp, and olefile.
This does not install services or write system directories. Signature verification
must precede extraction; see runtime/SOURCE.md for the official source.
"""
from pathlib import Path
import io
import shutil
import struct
import subprocess
import tempfile
import olefile

root = Path(__file__).resolve().parent
installer = root / '.tmp/msmpisetup.exe'
payload = installer.read_bytes()
destination = root / 'runtime/msmpi'
destination.mkdir(parents=True, exist_ok=True)

signature = subprocess.run([
    'powershell.exe', '-NoProfile', '-Command',
    f"(Get-AuthenticodeSignature -LiteralPath '{str(installer).replace(chr(39), chr(39)*2)}').Status"
], check=True, capture_output=True, text=True)
if signature.stdout.strip() != 'Valid':
    raise RuntimeError('Official installer signature is not Valid')

start = 0
for index in range(10):
    offset = payload.find(bytes.fromhex('d0cf11e0a1b11ae1'), start)
    if offset < 0:
        break
    start = offset + 8
    with olefile.OleFileIO(io.BytesIO(payload[offset:])) as archive:
        for stream in archive.listdir():
            if archive.get_size(stream) < 100000:
                continue
            data = archive.openstream(stream).read()
            if not data.startswith(b'MSCF'):
                continue
            # List each cabinet first. Only extract the x64 payload containing msmpi64.dll.
            cabinet = root / f'.tmp/mpi_payload_{index}.cab'
            cabinet.write_bytes(data)
            listing = subprocess.run(['tar', '-tf', str(cabinet)], check=True,
                                     capture_output=True, text=True).stdout
            if 'msmpi64.dll' not in listing:
                continue
            with tempfile.TemporaryDirectory(dir=root / '.tmp') as staging:
                subprocess.run(['tar', '-xf', str(cabinet), '-C', staging], check=True)
                for file in Path(staging).iterdir():
                    # Retain 32-bit copies separately; use 64-bit DLLs for Python.
                    name = {'msmpi.dll': 'msmpi_x86.dll',
                            'msmpires.dll': 'msmpires_x86.dll'}.get(file.name, file.name)
                    shutil.copy2(file, destination / name)
                for name in ['msmpi', 'msmpires']:
                    shutil.copy2(destination / f'{name}64.dll', destination / f'{name}.dll')
            print('Extracted x64 Microsoft MPI runtime to', destination)
            raise SystemExit(0)
raise RuntimeError('No x64 Microsoft MPI cabinet found')
