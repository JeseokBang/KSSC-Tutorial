function Area=ShoeLaceFormula(x1,y1,x2,y2,x3,y3)
    Area=1/2*abs(x1*y2+x2*y3+x3*y1-x2*y1-x3*y2-x1*y3);
end