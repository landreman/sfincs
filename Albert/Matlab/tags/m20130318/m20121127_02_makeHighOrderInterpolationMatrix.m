function m=m20121127_02_makeHighOrderInterpolationMatrix(x,y,L,GOrH)
% This function produces a matrix for interpolating a function from one
% grid (x) to another (y).  The interpolation is performed using cubic 
% functions.  This routine could be applied in many ways, but the 
% particular application it was designed for is interpolating Rosenbluth 
% potentials.  A variety of options are therefore available for 
% extrapolating the potentials off the right end of the x grid, depending 
% on which of the two potentials we are considering.  This function is not
% designed forextrapolating off the left end of the x grid, so use with 
% caution.
%
% Inputs:
% x = grid on which we know a function.
% y = new grid onto which we want to interpolate the function.
% L = Legendre mode number of the potential we are considering.  This
%     parameter is ignored when 'GOrH' = 'f'.
% GOrH = method for extrapolating to the right of the x grid.
%   Possible options for GOrH:
%   'f' = return 0.
%   'H' = extrapolate as 1/(y ^ (L+1))
%         (suitable for the 1st Rosenbluth potential.)
%   'G' = extrapolate using equation (27) of 
%         http://arxiv.org/pdf/1210.5289v1.pdf
%         (suitable for d^2/dv^2 of the 2nd Rosenbluth potential.)
%   Note: for the pure plasma code, it is not necessary to extrapolate the
%   potentials off the grid, so the 'L' and 'GOrH' parameters do not
%   matter.  They do matter in the impure-plasma code however.
%
% Outputs:
% m = interpolation matrix

m = zeros(numel(y), numel(x));
for i=1:numel(y)
    index = find(abs(x-y(i))<1e-12,1);
    if ~isempty(index)
        m(i,index)=1;
        
    else
        
        index = find(x>=y(i),1);
        if numel(index)>1
            error('numel(index)>1')
        end
        if numel(index)>0
            if (index>(numel(x)-2))
                indicesToUse = ((-3):0) + numel(x);
            elseif index<3
                indicesToUse = 1:4;
            else
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
            elseif strcmp(GOrH,'f')
                % Extrapolated value should be zero, so do nothing.
            else
                error('Invalid GOrH')
            end
        end
    end
end

end