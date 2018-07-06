function [mssd] = compute_mssd( x )

% this function computes the meas successive squared difference (MSSD).
% MSSD measures moment to moment changes in signal where the mean may be
% shifting.

% from http://onlinelibrary.wiley.com/store/10.1002/0471667196.ess2635.pub2/asset/ess2635.pdf?v=1&t=jcwju37j&s=76c5fff6ace51dde0f4d4f12c2aecfbb04a37c0a&systemMessage=Please+be+advised+that+we+experienced+an+unexpected+issue+that+occurred+on+Saturday+and+Sunday+January+20th+and+21st+that+caused+the+site+to+be+down+for+an+extended+period+of+time+and+affected+the+ability+of+users+to+access+content+on+Wiley+Online+Library.+This+issue+has+now+been+fully+resolved.++We+apologize+for+any+inconvenience+this+may+have+caused+and+are+working+to+ensure+that+we+can+alert+you+immediately+of+any+unplanned+periods+of+downtime+or+disruption+in+the+future.
% Successive Differences 

T = size(x,1);

% mssd
x_d = diff(x);
mssd = sum(x_d.^2);
mssd = mssd/(T-1);