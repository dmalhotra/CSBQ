function Lambda = Lambda_SBT(pan, logterm)
% LAMBDA_SBT   fill Nystrom discretization of "local" Lambda tensor in SBT
%
% L = Lambda_SBT(pan,logterm) returns 3N*3N matrix (where N is number of nodes
%  in the panelization pan), zero apart from the j'th 3*3 diagonal block being
%  the tensor
%                     (I-3shat.shat^T) + 2*(I+shat.shat^T).logterm
%
%  where shat = \hat{s} = unit tangent at j'th node, and user must supply
%  the number logterm = log(4*L/(pi*eps)) where L is perimeter of curve.
%  This is the classical SBT local term which hits force density values on
%  nodes, in the "xyz-fast, node slow" ordering.
%
%  The only test so far is DRAG_TORUS_EDGEWISE_SBT

% Barnett 1/18/22

tx = horzcat(pan.tx);                  % s-hat's at nodes, 3*N
N = size(tx,2);                        % # total nodes

% SBT "local" tensor: (I-3shat.shat^T) + 2*(I+shat.shat^T).logterm, not const
Lambda = (1+2*logterm)*eye(3*N);       % start w/ isotropic part of Lambda
for j=1:N, j3=3*j+(-2:0);              % loop over node inds, dyadic part
  shat = tx(:,j);                      % this node hat{s}, unit tangent
  Lambda(j3,j3) = Lambda(j3,j3) + (-3+2*logterm)*(shat*shat');  % 3x3 diag blk
end
