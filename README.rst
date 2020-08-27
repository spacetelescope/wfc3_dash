Introduction
------------

This package is intended to aid users in properly reducing their DASH data. While there are many portions of the pipeline that will need to be specialized depending on independent data/observations the general outline of steps should be useful to all users in their reduction. Further, the notebooks which provide a pipeline walk through will be generally applicable to all DASH data and provide a helpful interpretation of the best practices as described in Momcheva et. al (INSERT LINK).

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

  python setup.py [install|develop]


Lacosmicx Install
~~~~~~~~~~~~~~~~~
The last install needed is a package known as LACosmicx which is used within the pipeline to aid in cosmic ray rejection. The package can be downloaded and installed from the command line as follows: 

::

  git clone https://github.com/cmccully/lacosmicx.git
  cd lacosmicx
  python setup.py install
