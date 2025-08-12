function [button1, button2, button3] = gui_createbuttons(parentfig, in_sc_openscedlg)

    fig_pos = get(parentfig, 'Position'); 
    fig_width = fig_pos(3);
    fig_height = fig_pos(4);
    btn_width = 100; 
    btn_height = 25;
    btn_x = (fig_width - btn_width) / 2;
    btn_y = (fig_height - btn_height) / 1.618;

if gui.i_isuifig(parentfig)
    mfolder = fileparts(mfilename('fullpath'));
    iconfile = fullfile(mfolder, '..','assets', 'Images', 'box3d-center.jpg');
    
    button1 = uibutton(parentfig,"Text",'Import Data...', ...
        'ButtonPushedFcn', in_sc_openscedlg,...
        'Position', [btn_x btn_y btn_width+5 btn_height],...
        'Icon', iconfile);
    button2 = uilabel(parentfig, 'FontSize', 12,...
        'position', [btn_x btn_y+25 btn_width+50 btn_height],...
        'Text','Ready to explore.');
    if nargout==3
        button3 = uibutton(parentfig,"Text",'Switch to New SCGEATOOL', ...
            'ButtonPushedFcn', @in_askswitch,...
            'Position', [0 0 btn_width*2 btn_height]);
    end

else
    
    button1 = uicontrol('Parent', parentfig,...
        'Style', 'pushbutton',...
        'Units', 'pixels',...
        'Position', [btn_x btn_y btn_width btn_height],...
        'String', 'Import Data...',...
        'Callback', in_sc_openscedlg,...
        'ButtonDownFcn', in_sc_openscedlg,...
        'KeyPressFcn', in_sc_openscedlg, 'Tooltip','Click or Press i');
    button2 = uicontrol('Parent', parentfig,...
        'style','text',...
        'Units', 'pixels',...
        'position', [btn_x btn_y+25 btn_width btn_height],...
        'FontSize', 9,...    
        'string','Ready to explore.');
    if nargout==3
        button3 = uicontrol('Parent', parentfig,...
            'Style', 'pushbutton',...
            'Units', 'pixels',...
            'Position', [0 0 btn_width*1.5 btn_height],...
            'String', 'Switch to New SCGEATOOL',...
            'Callback', @in_askswitch,...
            'ButtonDownFcn', @in_askswitch);
    end
    if get(parentfig,'AutoResizeChildren')==0
        set(parentfig,'resizefcn',{@gui_myresizefun, button1, button2});
    end

   
end

    function in_askswitch(~, ~)
        group = "scgeatoolbox";
        pref = "switch2uifigure";
        % strcmp(getpref(group,pref), 'ask')
        setpref(group, pref, 'ask');
        gui.myHelpdlg(parentfig,'You will be prompted to switch the next time you launch SCGEATOOL.');
    end

end


function gui_myresizefun(src, ~, butt1, butt2)
    fig_pos = get(src, 'Position');
    fig_width = fig_pos(3);
    fig_height = fig_pos(4);

    btn_width = 100;
    btn_height = 25;
    btn_x = (fig_width - btn_width) / 2;
    btn_y = (fig_height - btn_height) / 1.618;
    
    set(butt1, 'Position', [btn_x btn_y btn_width btn_height]);
    set(butt2, 'Position', [btn_x btn_y+25 btn_width btn_height]);
end

