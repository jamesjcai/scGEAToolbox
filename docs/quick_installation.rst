Quick installation
==================

Run the following code in `MATLAB`:

.. code-block:: matlab

  tic
  disp('Installing scGEAToolbox...')
  unzip('https://github.com/jamesjcai/scGEAToolbox/archive/master.zip');
  addpath('./scGEAToolbox-master');  
  toc
  if exist('cdgea.m','file')
      disp('scGEAToolbox installed!')
  end
  savepath;
  
