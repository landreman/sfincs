function m = m20160925_01_periodicInterpolation(x,y,period,stencil)
% x = old grid
% y = new grid

% Stencil can be 1, 2, 3, 4, or 5.

% x should be monotonically increasing, and in the range [0,period):
assert(min(x) >= 0)
assert(max(x) < period)
for j = 2:numel(x)
    assert(x(j-1) < x(j))
end

y = mod(y,period);

m = zeros(numel(y), numel(x));
N = numel(x);

switch stencil
    case -1
        % Nearest-neighbor
        for i=1:numel(y)
            [residual1, index1] = min(abs(x-y(i)));
            [residual2, index2] = min(abs(x+period-y(i)));
            [residual3, index3] = min(abs(x-period-y(i)));
            [~,index4] = min([residual1,residual2,residual3]);
            indices = [index1,index2,index3];
            m(i,indices(index4))=1;
        end
        
    case 1
        % Nearest-neighbor
        xc = x(:);
        x3 = [xc-period; xc; xc+period];
        for i=1:numel(y)
            [~, index] = min(abs(x3-y(i)));
            m(i,mod(index-1,N)+1)=1;
        end
        
    case 2
        for i=1:numel(y)
            index = find(x>=y(i),1); % Find just the 1st index
            if numel(index)>1
                index
                error('Should not get here: numel(index)>1')
            end
            
            x0=0;
            x1=0;
            
            if numel(index)==0
                % y(i) is to the right of all the x's
                indicesToUse = [N, 1];
                x1 = period;
            elseif index==1
                % y(i) is to the left of all the x's
                indicesToUse = [N, 1];
                x0 = -period;
            else
                % y(i) is inside the range of x
                indicesToUse = [-1,0] + index;
            end
            
            x0 = x0 + x(indicesToUse(1));
            x1 = x1 + x(indicesToUse(2));
            xx = y(i);
            
            m(i,indicesToUse(1)) = (xx-x1)/( (x0-x1));
            m(i,indicesToUse(2)) = (xx-x0)/( (x1-x0));
            
        end
        
    case 3
        % First find the nearest neighbor
        xc = x(:);
        x3 = [xc-period; xc; xc+period];
        for i=1:numel(y)
            [~, index] = min(abs(x3-y(i)));
            centerIndex = mod(index-1,N)+1;
            indicesToUse = [-1,0,1]+centerIndex;
            
            x0=0;
            x1=0;
            x2=0;
            if centerIndex==1
                indicesToUse(1) = N;
                x0 = -period;
            elseif centerIndex==N
                indicesToUse(3) = 1;
                x2 = period;
            end

            x0 = x0 + x(indicesToUse(1));
            x1 = x1 + x(indicesToUse(2));
            x2 = x2 + x(indicesToUse(3));
            xx = y(i);
            
            m(i,indicesToUse(1)) = (xx-x1)*(xx-x2)/( (x0-x1)*(x0-x2));
            m(i,indicesToUse(2)) = (xx-x0)*(xx-x2)/( (x1-x0)*(x1-x2));
            m(i,indicesToUse(3)) = (xx-x0)*(xx-x1)/( (x2-x0)*(x2-x1));
        end


    case 4
        for i=1:numel(y)
            index = find(x>=y(i),1); % Find just the 1st index
            if numel(index)>1
                index
                error('Should not get here: numel(index)>1')
            end
            
            x0=0;
            x1=0;
            x2=0;
            x3=0;
            
            if numel(index)==0
                % y(i) is to the right of all the x's
                indicesToUse = [N-1, N, 1, 2];
                x2 = period;
                x3 = period;
            elseif index == N
                indicesToUse = [N-2, N-1, N, 1];
                x3 = period;
            elseif index==1
                % y(i) is to the left of all the x's
                indicesToUse = [N-1, N, 1, 2];
                x0 = -period;
                x1 = -period;
            elseif index==2
                indicesToUse = [N, 1, 2, 3];
                x0 = -period;
            else
                indicesToUse = ((-2):1) + index;
            end
            
            x0 = x0 + x(indicesToUse(1));
            x1 = x1 + x(indicesToUse(2));
            x2 = x2 + x(indicesToUse(3));
            x3 = x3 + x(indicesToUse(4));
            xx = y(i);
            
            m(i,indicesToUse(1)) = (xx-x1)*(xx-x2)*(xx-x3)/( (x0-x1)*(x0-x2)*(x0-x3));
            m(i,indicesToUse(2)) = (xx-x0)*(xx-x2)*(xx-x3)/( (x1-x0)*(x1-x2)*(x1-x3));
            m(i,indicesToUse(3)) = (xx-x0)*(xx-x1)*(xx-x3)/( (x2-x0)*(x2-x1)*(x2-x3));
            m(i,indicesToUse(4)) = (xx-x0)*(xx-x1)*(xx-x2)/( (x3-x0)*(x3-x1)*(x3-x2));
            
        end
        
    otherwise
        error('Invalid stencil: %d',stencil)
end