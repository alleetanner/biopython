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

## Anaconda and Jupyter Lab
To provide a consistent interface, we will be using Jupyter Lab as our workspace. If you are familiar with other development platforms, for example [VSCode](https://code.visualstudio.com/) or [PyCharm](https://www.jetbrains.com/pycharm/), you are welcome to use these, but be sure you can also handle the other requirements for your set up (ie, how to install packages and how to access the terminal). 

If you haven't already, download and install [Anaconda](https://www.anaconda.com/). Anaconda is a data science platform for working with Python and R (as well as other things). It includes Jupyter Lab (a browser-based environment providing _development tools_ such as a text editor, terminal, and Python console), and Jupyter Notebooks (a browser-based IDE for creating "Notebooks": _hybrid documents of runnable code, with Markdown_ (ie, commentary) blocks to document what the code is doing and why). We will be working in Jupyter Lab - from the Anaconda Navigator starting screen, click the Jupyter Lab button and it will start a new tab in your default browser.

For this workshop we will need both a text editor and a terminal pane in this tab. You can open these by clicking on the associated icons on the start screen, or by going `File` > `New` > `Text file` and `File` > `New` > `Terminal`. Arrange these as you like by dragging the tabs - we usually teach with a text editor on the left and a terminal on the right. Also, create a new folder to work in, by navigating using the folder icon in the menu on the left and creating one. Make sure you move into this folder in your terminal, and when you save your text file ensure it is in this folder.

## Installing BioPython
Finally, we need to download BioPython and install it as part of your local Python. We will do this both for Anaconda, and in the terminal - they may be using different installs of Python, so we will install it both ways to be sure!

### BioPython for Anaconda
In your Anaconda Navigator window (the one with the launch buttons, not the tab in your browser), click the `Environments` tab to the left. You should have `base(root)` selected in the second column. On the right column, from the drop-down select `Not installed`, and type `biopython` in the search box to the right. Select the BioPython package, click `Apply`, and then confirm the install. This might take a couple of minutes.

### BioPython in the terminal
If you are using Windows, the terminal that Jupyter provides should be running the install of Python that we just installed BioPython for. However sometimes (especially on MacOS), it might not be using the same Python install. If you are having issues running BioPython in the Jupyter terminal, try installing BioPython there: go to the terminal and enter

```shell
pip install biopython
```

and follow the prompts to download and install the package. We should be ready to go now!

{{< admonition type="info" title="Package managers" open=true >}}
`pip` is a package manager for Python. `conda` is the package manager for Anaconda. A package manager is a program which retrieves extra Python libraries, installs, updates, and generally organises your Python setup. Whenever we run an `import` command, Python will look for modules and libraries, first in your working folder (where you might have saved your own modules), then in Python's local folders of installed packages (ie, back near your operating system file - so you won't see them unless you go looking for them!). Apart from Python standard libraries (like `math` and `sys` for example), we have to download and install anything else (sometimes called _3rd party libraries_: things not made by you, nor included as base Python, but made by other developers).

That is where `pip` (or `conda`) comes it. It will search package repositories, typically [PyPI](https://pypi.org/), download the appropriate package, and install it, along with any other dependencies it might need. Here are some further `pip` commands, to examine what you have currently installed, etc:
```shell
pip list
pip info
pip help
```
{{< /admonition >}}
