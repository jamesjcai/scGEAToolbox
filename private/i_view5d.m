function i_view5d(s5d)
    if size(s5d,2)<5
        s5d=[s5d,zeros(size(s5d,1),5-size(s5d,2))];
    end
    i_view5dIDX=s5d(:,3);
    assignin('base','i_view5dIDX',i_view5dIDX);
    subplot(1,2,1);
    h1=scatter3(s5d(:,1),s5d(:,2),s5d(:,3),10);
    xlabel('dim 1')
    ylabel('dim 2')
    zlabel('dim 3')
    grid on
    box on
    h1.ZDataSource='i_view5dIDX';

    subplot(1,2,2);
    h2=scatter3(s5d(:,3),s5d(:,4),s5d(:,5),10);
    xlabel('dim 3')
    ylabel('dim 4')
    zlabel('dim 5')
    grid on
    box on
    h2.XDataSource='i_view5dIDX';
    hLD = linkdata('on');
    evalin('base','h=findobj(gcf,''type'',''axes'');');
    evalin('base','hlink = linkprop(h,{''CameraPosition'',''CameraUpVector''});');
    evalin('base','rotate3d on');
end
