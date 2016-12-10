function m=m20120925_09_makeHighOrderUniformRegriddingMatrix(x,y,L,GOrH)
% x = old grid
% y = new grid

m = zeros(numel(y), numel(x));
for i=1:numel(y)
    index = find(abs(x-y(i))<1e-12,1);
    if ~isempty(index)
        m(i,index)=1;
        %fprintf('Adding a 1 at row %d, column %d\n',i, index)
    else
        
        index = find(x>=y(i),1);
        if numel(index)>1
            index
            error('numel(index)>1')
        end
        if numel(index)>0
            if (index>(numel(x)-2))
                %fprintf('A\n')
                indicesToUse = ((-3):0) + numel(x);
            elseif index<3
                indicesToUse = 1:4;
            else
                %fprintf('B\n')
                indicesToUse = ((-2):1) + index;
            end
            
            x0 = x(indicesToUse(1));
            x1 = x(indicesToUse(2));
            x2 = x(indicesToUse(3));
            x3 = x(indicesToUse(4));
            xx = y(i);
            
            m(i,indicesToUse(1)) = (xx-x1)*(xx-x2)*(xx-x3)/( (x0-x1)*(x0-x2)*(x0-x3));
            m(i,indicesToUse(2)) = (xx-x0)*(xx-x2)*(xx-x3)/( (x1-x0)*(x1-x2)*(x1-x3));
            m(i,indicesToUse(3)) = (xx-x0)*(xx-x1)*(xx-x3)/( (x2-x0)*(x2-x1)*(x2-x3));
            m(i,indicesToUse(4)) = (xx-x0)*(xx-x1)*(xx-x2)/( (x3-x0)*(x3-x1)*(x3-x2));
        else
            % We're extrapolating.
            if strcmp(GOrH,'H')
                % H and x dH/dx behave as D \propto 1/(x ^ (L+1)).
                m(i,end) = (x(end)/y(i))^(L+1);
            elseif strcmp(GOrH,'G1')
                % Consider only 1 of the 2 polynomial terms in G:
                % G \propto 1/x^(L-1)
                m(i,end) = (x(end)/y(i))^(L-1);
            elseif strcmp(GOrH,'G')
                xi = x(end);
                xj = x(end-1);
                denomQ = xi^2-xj^2;
                denomP = 1/xi^2 - 1/xj^2;
                m(i,end) = xi^(L+1)/denomP/y(i)^(L+3)  + xi^(L+3)/denomQ/y(i)^(L+1);
                m(i,end-1) = -xj^(L+1)/denomP/y(i)^(L+3)  - xj^(L+3)/denomQ/y(i)^(L+1);
                
                % Note: the lines below are wrong:
                %m(i,end) = xi^(L-1)/denomP/y(i)^(L+1)  + xi^(L+1)/denomQ/y(i)^(L-1);
                %m(i,end-1) = -xj^(L-1)/denomP/y(i)^(L+1)  - xj^(L+1)/denomQ/y(i)^(L-1);
            elseif strcmp(GOrH,'f')
                % Extrapolated value should be zero, so do nothing.
            else
                error('Invalid GOrH')
            end
        end
    end
end

end