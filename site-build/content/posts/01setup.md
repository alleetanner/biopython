---
title: "• Set up and requirements"
subtitle: "Creating a suitable environment for running BioPython"

date: 2022-11-01T00:00:00+01:00

fontawesome: true
linkToMarkdown: true

toc:
  enable: true
  keepStatic: false
  auto: false
code:
  copy: true
math:
  enable: true
share:
  enable: false
comment:
  enable: false
---

### Anaconda and Jupyter Lab
In order to provide a consistent interface, we will be using Jupyter Lab as our workspace. If you are familiar with other tools, such as PyCharm or VSCode, you are welcome to use these, but be sure you can also handle the specific requirements for your set up (ie, how to install packages and how to access the terminal). 

If you haven't already, download [Anaconda](https://www.anaconda.com/). Anaconda is a platform providing data science tools for Python and R. The standard install includes things like Jupyter Notebooks (a browser-based IDE for creating hybrid documents of runnable code, with Markdown blocks to document what the code is doing), and Jupyter Lab (a browser-based environment providing a text editor, terminal, and Python console). We will be working in Jupyter Lab - from the Anaconda Navigator starting screen, click the Jupyter Lab button and it will start a session in a new tab in your default browser.

![The Anaconda launch window](/assets/jupyter-lab-screen.png#center)

For this workshop we will need both a text editor and a terminal pane in this tab. You can open these by clicking on the associated icons on the start screen, or by going `File` > `New` > `Text file` and `File` > `New` > `Terminal`. Arrange these as you like by dragging the tabs - we usually teach with a text editor on the left and a terminal on the right. Also, created a new folder to work in, by navigating using the folder icon in the menu on the left and creating one. Make sure you move into this folder in your terminal, and when you save your text file ensure it is in this folder.


### Installing BioPython
Finally, we need to download BioPython and install it as part of your local Python. We will do this both for Anaconda, and in the terminal - they may be using different installs of Python, so we will install it both ways to be sure!

#### BioPython for Anaconda
In your Anaconda Navigator window (the one with the launch buttons, not the tab in your browser), click the `Environments` tab to the left. You should have `base(root)` selected in the second column. On the right column, from the drop-down select `Not installed`, and type `biopython` in the search box to the right. Select the BioPython package, click `Apply`, and then confirm the install. This might take a couple of minutes.

![The Anaconda Environment window]([https://github.com/alleetanner/biopython/blob/main/assets/anaconda-env-screen.png#center)

#### BioPython in the terminal
Admittedly, this can be confusing, but your terminal might be using its own version of Python. So, we should install BioPython there too. In your Jupyter Lab tab in your browser, go to the terminal and enter
```
pip install biopython
```
and follow the prompts to download and install the package. We should be ready to go now!

{{< admonition type="info" open=true >}}
`pip` is a package manager for Python. Whenever we run an `import` command in Python, it will look for modules and libraries to import, first in your working folder (where you might have put your own modules), then in Python's local repository of installed packages. Apart from standard packages fundamental to Python, we have to download and install anything else (sometimes called 3rd party packages).

That is where `pip` comes it. It will search package repositories, typically [PyPI](https://pypi.org/), download the appropriate package, and install it, along with any other dependencies it might need.
{{< /admonition >}}
