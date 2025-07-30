Install as a MATLAB Add-On (Recommended)
=========================================

1. **Open the Add-On Explorer**
   - In MATLAB, go to the **Home** tab.
   - Click on the **Add-Ons** icon to open the Add-On Explorer.

2. **Search for scGEAToolbox**
   - In the search bar, type **"scGEAToolbox"** and press Enter.

3. **Select and Install**
   - Find **"scGEAToolbox (single-cell Gene Expression Analysis Toolbox)"** in the results.
   - Click the **Add** button to install the toolbox.

4. **Launch scGEAToolbox**
   - To start using scGEAToolbox, enter:

     .. code-block:: matlab

        scgeatool

Quick Installation (For Developers)
===================================

Run the following code in MATLAB:

.. code-block:: matlab

   unzip('https://github.com/jamesjcai/scGEAToolbox/archive/main.zip');
   addpath('./scGEAToolbox-main');
   savepath(fullfile(userpath,'pathdef.m'));
   % savepath;