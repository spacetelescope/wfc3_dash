THIS VERSION IS NOW DEPRECATED - TO FIND THE CURRENT DASH NOTEBOOK PLEASE GO TO: 
---------------------------------------------------------------------------------

https://github.com/spacetelescope/WFC3Library
You will need to install the `wfc3_env` with conda and the `lacosmicx` as well using the instructions available on the main WFC3Library page and the Dash notebook `read me`. 






Introduction
------------

This package is intended to aid users in properly reducing their DASH data. While there are many portions of the pipeline that will need to be specialized depending on independent data/observations the general outline of steps should be useful to all users in their reduction. Further, the notebooks which provide a pipeline walk through will be generally applicable to all DASH data and provide a helpful interpretation of the best practices as described in `Momcheva et. al 2016 <https://arxiv.org/pdf/1603.00465.pdf>`_. 

An in-depth Instrument Science Report that walks through this package and the accompanying Jupyter Notebook tutorial is now available `here  <https://www.stsci.edu/files/live/sites/www/files/home/hst/instrumentation/wfc3/documentation/instrument-science-reports-isrs/_documents/2021/2021-02.pdf>`_.

Installation
------------

Download Anaconda or Miniconda
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You must first have a working installation of ``anaconda`` or ``miniconda`` for Python 3.  If you do not yet have this on your system, you can visit the following links for download and installation instructions:

- `Anaconda <https://www.anaconda.com/download/>`_
- `Miniconda <https://conda.io/en/latest/miniconda.html>`_

Obtain the ``wfc3_dash`` Package
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To obtain ``wfc3_dash`` clone the repository directly from GitHub:

::

  git clone https://github.com/spacetelescope/wfc3_dash.git
  cd wfc3_dash

Environment Installation
~~~~~~~~~~~~~~~~~~~~~~~~
First, to make sure your ``conda`` is up to date:

::

  conda update conda


Next, you should activate the ``base`` environment:

::

  conda activate base


Then you can create the ``dash`` ``conda`` environment with the following command:

::

  conda create -n dash stsci notebook astroquery numpy=1.17.4


where the numpy version is specifed to the current version that works with `Tweakreg`, a vital package to reducing DASH data. 

From there a user can activate the `dash` environment at any point with:

::

  conda activate dash


Package Installation
~~~~~~~~~~~~~~~~~~~~

In order to install the ``dash`` package within the ``conda`` environment, run the `dash` setup script when in the `wfc3_dash` main directory:

::

  python setup.py install
 

Package Updates
~~~~~~~~~~~~~~~

There will be periodic changes to the ``dash`` package as we update and expand the codes and notebook examples based on current plans and user input. To update your `dash` package, you will first go to the `wfc3_dash` main directory and check that you are on the `master` branch:

::

  git branch

If you are not on the master branch, you can switch back using: 

:: 

  git checkout master
  
From there you can pull any changes to the code or notebooks using:

:: 

  git pull

And install those changes using the setup script again: 

::

  python setup.py install


Lacosmicx Install
~~~~~~~~~~~~~~~~~
The last install needed is a package known as LACosmicx which is used within the pipeline to aid in cosmic ray rejection. The package can be downloaded and installed from the command line as follows: 

::

  git clone https://github.com/cmccully/lacosmicx.git
  cd lacosmicx
  python setup.py install
