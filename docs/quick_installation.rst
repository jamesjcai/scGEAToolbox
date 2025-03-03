Install via MATLAB Add-On (`.mltbx`) (Recommended)
==================================================
1. **Download** `scGEAToolbox.mltbx` from:

   - 🔗 [`GitHub Releases <https://github.com/jamesjcai/scGEAToolbox/releases>`__]

   - 🔗 [`MATLAB File Exchange <https://www.mathworks.com/matlabcentral/fileexchange/72917-scgeatoolbox-single-cell-gene-expression-analysis-toolbox>`__]

2. **Install** the toolbox:  
   - **Double-click** `scGEAToolbox.mltbx`  
   - MATLAB will open the **Add-On Manager**  
   - Click **Install**  

3. **Verify the installation** by running:  
  
.. code-block:: matlab
  
  matlab.addons.installedAddons
  
Ensure scGEAToolbox appears in the list.


Quick Installation (For Developers)
==================================

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
  
