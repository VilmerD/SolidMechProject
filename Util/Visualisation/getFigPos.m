function figconfig = getFigPos()
wboxes = 2;     hboxes = 2;
wfiggrid = 2;   hfiggrid = 1;

resolution = get(0, 'screensize');
w = resolution(3);     h = resolution(4);
boxwidth = w/wboxes;     boxheight = h/hboxes;

padding = 0.025;
wpadding = ceil(padding*w);   hpadding = ceil(padding*h);

figwidth = (w - wpadding*(1 + wfiggrid*wboxes))/(wfiggrid*wboxes);
figheight = (h - hpadding*(1 + hfiggrid*hboxes))/(hfiggrid*hboxes);
figwidth = floor(figwidth);     figheight = floor(figheight);


figconfig = cell(hboxes, wboxes);
x0 = wpadding; 
y0 = h - (hpadding + figheight);
for k = 1:wboxes
    xindex = k;
    x = x0 + (k - 1) * (boxwidth - wpadding/2);
    for i = 1:hboxes
        yindex = i;
        figgridik = zeros(wfiggrid*hfiggrid, 4);
        y = y0 - (i - 1) * (boxheight + hpadding/2);
        for j = 1:wfiggrid
            xj = x + (j - 1)*(figwidth + wpadding);
            for l = 1:hfiggrid
                yj = y + (l - 1)*(figheight + hpadding);
                
               figgridik(j + l - 1, :) = [xj, yj, figwidth, figheight];
            end
        end
        figconfig{xindex, yindex} = figgridik;
    end
end
end

