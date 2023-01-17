function plotLABpixel(X,averagePixelDepth)

n = length(X)-1;
hold on
for I = 1:n
    xx = [X(I) X(I+1)];
    yy = averagePixelDepth(I) * [1 1];
    plot(xx,yy,'Color',	'#7E2F8E','LineStyle','-.','LineWidth',3);
end

end