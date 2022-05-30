import argparse
import subprocess
import os


def main():
    parser = argparse.ArgumentParser(
        description=('converts .ipynb files to HTML with nbconvert and modifies'
                     ' for my website')
    )
    parser.add_argument(
        'file',
        help='the .ipynb file to convert to HTML'
    )
    args = parser.parse_args()

    if args.file.split('.')[1] != 'ipynb':
        raise ValueError('file must be a .ipynb file')
    
    proc = subprocess.run(['jupyter', 'nbconvert', '--to', 'html', '--HTMLExporter.theme=dark', args.file],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.DEVNULL)
    if proc.returncode != 0:
        raise ValueError('jupyter nbconvert failed')

    new_html = f'{args.file.split(".")[0]}.html'
    old_html = f'{args.file.split(".")[0]}_old.html'
    os.rename(new_html, old_html)

    bootstrap_trigger = '<meta name="viewport" content="width=device-width, initial-scale=1.0">'
    bootstrap_addition = '<link href="https://cdn.jsdelivr.net/npm/bootstrap@5.2.0-beta1/dist/css/bootstrap.min.css" rel="stylesheet" integrity="sha384-0evHe/X+R7YkIZDRvuzKMRqM+OrBnVFBL6DOitfPri4tjfHxaWutUpFmBp4vmVor" crossorigin="anonymous">\n'
    header_scroll_trigger = '<body class="jp-Notebook" data-jp-theme-light="false" data-jp-theme-name="JupyterLab Dark">'
    header_addition = """
    <div class="container text-white text-center">
      <main class="px-3">
        <h1>Molecule Club</h1>
      </main>
      <header class="d-flex flex-wrap justify-content-center py-3 mb-4 border-bottom">
        <ul class="nav nav-pills">
          <li class="nav-item px-1"><a href="/index.html" class="nav-link text-bg-light">Home</a></li>
          <li class="nav-item px-1"><a href="/resume.html" class="nav-link text-bg-light">Resume</a></li>
          <li class="nav-item px-1"><a href="/notebooks/index.html" class="nav-link text-bg-light">Notebooks</a></li>
        </ul>
      </header>
    </div>
    """
    scroll_replace = '<body class="jp-Notebook" data-jp-theme-light="false" data-jp-theme-name="JupyterLab Dark" style="overflow:auto;">'

    with open(old_html) as input_file:
        with open(new_html, 'w') as output_file:
            for line in input_file.readlines():
                if bootstrap_trigger in line:
                    output_file.write(bootstrap_addition)
                if header_scroll_trigger in line:
                    output_file.write(header_addition)
                    line = scroll_replace
                output_file.write(line)

if __name__ == '__main__':
    main()