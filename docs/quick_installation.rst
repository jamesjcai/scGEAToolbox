Quick installation
==================

Run the following code in `MATLAB`:

.. code-block:: matlab

  tic
  disp('Installing scGEAToolbox...')
  unzip('https://github.com/jamesjcai/scGEAToolbox/archive/main.zip');
  addpath('./scGEAToolbox-main');  
  toc
  if exist('scgeatool.m','file')
      disp('scGEAToolbox installed!')
  end  
  savepath(fullfile(userpath,'pathdef.m')); 
  % savepath;
  
