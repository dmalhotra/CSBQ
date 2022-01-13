function b = pan_brkpts(W)
% PAN_BRKPTS  sum up quadrature weights in panels and cumsum to breakpoints
%
% b = pan_brkpts(W) where W is a npan*p matrix with row j the quadrature
%  weights for panel j, and returns a column b of breakpoints with respect to
%  whatever parameter (needn't be arclength) the quadrature is for. b(1)=0
%  and b(j+1)-b(j) is the sum of the jth row of W. b(npan+1) is the total
%  parameter length.
%
% Called with no arguments does self-test.

% Barnett 1/12/22
if nargin==0, test_pan_brkpts; return; end
b = [0; cumsum(sum(W,2))];

%%%%%%
function test_pan_brkpts   % test round-trip from tpan -> quad wei -> tpan
p = 12;                    % order
npan = 16;
L = 2*pi;                  % period
tpan = L*(0:npan)'/npan;   % pan param breakpoints (first=0, last=2pi)
tpan(2:npan) = tpan(2:npan) + 0.2*L*(2*rand(npan-1,1)-1)/npan;  % unequal panels
pan = setup_pans(tpan,p);
W = horzcat(pan.v)';       % setup input mat (v are col vecs)
b = pan_brkpts(W);
norm(b-tpan)
