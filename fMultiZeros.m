function ListZeros = fMultiZeros(f,xmin,xmax,dx)
% ListZeros = fMultiZeros(f,xmin,xmax,dx) returns the list ListZeros of zeros of 
% function f between xmin and xmax. ListZeros is the empty list if no zero
% is found.


    ListZeros = [];
    x2=xmin;
    y2=f(x2);

    while x2<xmax
        x1=x2;
        y1=y2;
        x2=x2+dx;
        y2=f(x2);
        if (y1*y2<=0)                              
            ListZeros=[ListZeros,fzero(f,[x1,x2])]; 
        end
    end


