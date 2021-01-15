.. _getting_started:

Getting Started
===============

Qick installation
-----------------
Run the following code in `MATLAB`:

.. codeblock::

  tic
  disp('Installing scGEAToolbox...')
  unzip('https://github.com/jamesjcai/scGEAToolbox/archive/master.zip');
  addpath('./scGEAToolbox-master');
  toc
  if exist('cdgea.m','file')
      disp('scGEAToolbox installed!')
  end

Run DEMO_GETTING_STARTED
------------------------

.. codeblock::

 demo_Getting_Started


Using SC_SCATTER
----------------
For a quick exploratry data analysis using `SC_SCATTER` function

.. codeblock::

  cdgea;
  load example_data\testXgs.mat
  sc_scatter(X,s,g)

If everything goes right, you will see the figure like this:


