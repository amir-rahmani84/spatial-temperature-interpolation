"""Compile the authoritative LaTeX report; never rewrite README or report source.

Run from any directory: python prepare_documents.py
Requires pdfLaTeX (MiKTeX or TeX Live), separate from Python dependencies.
"""
from pathlib import Path
import argparse
import re
import shutil
import subprocess

ROOT=Path(__file__).resolve().parent

def compile_report(output_dir=None):
    executable=shutil.which('pdflatex')
    if not executable:
        raise FileNotFoundError('pdfLaTeX is required. Install MiKTeX or TeX Live and add pdflatex to PATH.')
    report=ROOT/'report'; source=report/'research_report.tex'
    if not source.is_file(): raise FileNotFoundError(source)
    output=Path(output_dir).resolve() if output_dir else report
    output.mkdir(parents=True,exist_ok=True)
    command=[executable,'-interaction=nonstopmode','-halt-on-error','-file-line-error',f'-output-directory={output}',source.name]
    for _ in range(2):
        completed=subprocess.run(command,cwd=report,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,encoding='utf-8',errors='replace')
        if completed.returncode:
            if completed.returncode in (3221225781,-1073741515):
                raise RuntimeError('pdfLaTeX cannot start: Windows error 0xC0000135 (missing runtime DLL). Repair the TeX installation or use a working pdfLaTeX installation; no report was regenerated.')
            raise RuntimeError('pdfLaTeX failed; inspect research_report.log in the output directory.\n'+completed.stdout[-5000:])
    log=(output/'research_report.log').read_text(encoding='utf-8',errors='replace')
    problems=re.findall(r'Overfull [^\n]+|LaTeX Warning: (?:Reference|Citation)[^\n]+|LaTeX Warning: There were undefined references[^\n]*',log)
    if problems:
        raise RuntimeError('Report requires layout/reference review:\n'+'\n'.join(problems))
    print(f'Compiled {output / "research_report.pdf"}; README and LaTeX source unchanged.')

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,help='Optional directory for the PDF and TeX build files')
    args=parser.parse_args(); compile_report(args.output_dir)
