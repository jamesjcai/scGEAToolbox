function callback_Quantis(src, ~)
    [FigureHandle, ~] = gui.gui_getfigsce(src);

listitems = natsort(["Cisplatin","Carboplatin","Paclitaxel","Docetaxel","Doxorubicin","Cyclophosphamide","Methotrexate","5-Fluorouracil (5-FU)","Capecitabine","Gemcitabine","Vincristine","Etoposide","Imatinib","Trastuzumab","Bevacizumab","Rituximab","Pemetrexed","Nivolumab","Pembrolizumab","Tamoxifen"]);
promptstr = "Select a drug:";

if gui.i_isuifig(FigureHandle)
        [indx, tf] = gui.myListdlg(FigureHandle, listitems, ...
            promptstr);
    else     
        [indx, tf] = listdlg('PromptString', ...
            {promptstr}, ...
            'SelectionMode', 'single', 'ListString', listitems, ...
            'ListSize', [220 300], 'Name', ' ');
    end    

    if tf == 1        
        gui.myHelpdlg(FigureHandle, "You selected "+listitems(indx)+".");
    end
end