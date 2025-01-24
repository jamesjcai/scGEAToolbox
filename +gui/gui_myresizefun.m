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
