function [ m, V, invV ] = dense_moments( x_vals, pdf )
  dx = x_vals(2) - x_vals(1);
  if numel(pdf)==numel(x_vals)
    m = x_vals * pdf' * dx;
    Exx = x_vals.^2 * pdf' * dx;
    V = Exx - m.^2;
    invV = 1/V;
  else
    [X,Y] = meshgrid(x_vals);    
    XY = [ X(:), Y(:) ];    
    m = sum( bsxfun(@times, XY, pdf(:)) ) * dx * dx;
    m = reshape(m,2,1);    
    XY_outer = bsxfun(@times, reshape(XY', 2,1,[] ), reshape(XY', 1,2,[]) );
    Exx = sum( bsxfun(@times, XY_outer, reshape(pdf(:),1,1,[])), 3 ) * dx * dx;     
    
%     Exx = zeros(2,2);
%     for i=1:size(XY,1)
%       this_x = XY(i,:);      
%       Exx = Exx + this_x' * this_x * pdf(i) * dx * dx;
%     end

    
    V = Exx - m*m';
    invV = inv(V);
  end
end