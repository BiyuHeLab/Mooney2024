function h = draw_circle(x,y,r)
    hold on
    % d = r*2;
    % px = x-r;
    % py = y-r;
    % h = rectangle('Position',[px py d d],'Curvature',[1,1]);
    % daspect([1,1,1])

    th = 0:pi/50:2*pi;
    xunit = r * cos(th) + x;
    yunit = r * sin(th) + y;
    h = plot(xunit, yunit);
    %from https://www.mathworks.com/matlabcentral/answers/98665-how-do-i-plot-a-circle-with-a-given-radius-and-center
    hold off