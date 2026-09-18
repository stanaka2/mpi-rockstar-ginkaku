# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'MPI-Rockstar'
copyright = '2024, Tomoyuki Tokuue, Tomoaki Ishiyama, Ken Osato, Satoshi Tanaka, and Peter Behroozi'
author = 'Tomoaki Ishiyama'
release = '1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'myst_parser',
#    'sphinxcontrib.mermaid',
#    'sphinx_markdown_tables',
#    'sphinx_copybutton',
#    'sphinx_diagrams',
]

templates_path = ['_templates']
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
#html_theme = 'furo'
html_static_path = ['_static']
#html_css_files = ['custom.css']
html_show_sourcelink = False

source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}


myst_enable_extensions = [
    "dollarmath",        
    "amsmath",           
    "colon_fence",       
    "deflist",           
    "html_admonition",   
    "html_image",        
    "linkify",
    "strikethrough",
]

myst_heading_anchors = 3
html_css_files = ['custom.css']
html_show_sourcelink = False

source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}


myst_enable_extensions = [
    "dollarmath",        
    "amsmath",           
    "colon_fence",       
    "deflist",           
    "html_admonition",   
    "html_image",        
    "linkify",
    "strikethrough",
]

myst_heading_anchors = 3
