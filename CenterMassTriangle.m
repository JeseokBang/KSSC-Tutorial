function [xg,yg,wg]=CenterMassTriangle(x1,y1,w1,x2,y2,w2,x3,y3,w3)
    xg=(x1+x2+x3)/3;
    yg=(y1+y2+y3)/3;
    wg=(w1+w2+w3)/3;
%     xg=(w1*x1+w2*x2+w3*x3)/(w1+w2+w3);
%     yg=(w1*y1+w2*y2+w3*y3)/(w1+w2+w3);
%     wg=(w1+w2+w3)/3;
end