function [ W ] = make_W_connectivity(isold, isupstream )
%MAKE_W_CONNECTIVITY creates W - the connectivity matrix for proprioception
%       Input:  isold = 1 for old W, =0 for new W
%               isupstream = 1 for upstream coupling, =0 for downstream

global eps alpha beta cee dim;

if isold == 1
    %%%--- OLD W
    %proprioceptive coupling matrix
    W = zeros(dim,dim);
    a_m = cee/(1+alpha);
    if isupstream==1
        % posterior->anterior coupling (upstream)
        b_up = alpha*a_m/(1+beta);
        b_down = beta*b_up;
    else
        % anterior->posterior coupling (downstream)
        b_down = alpha*a_m/(1+beta);
        b_up = beta*b_down;
    end
    W(1,2) = b_up;
    W(1,1) = cee-W(1,2);
    for ii=2:dim-1
        W(ii, ii) = a_m;
        W(ii,ii-1) = b_down;
        W(ii,ii+1) = b_up;
    end
    W(dim,dim-1) = b_down;
    W(dim,dim) = cee-W(dim,dim-1);
else
    %%%--- NEW W
    W = zeros(dim,dim);
    w_m = cee;
    if isupstream==1
        %posterior->anterior coupling (upstream)
        w_down = beta*eps;
        w_up = eps;
    else
        %anterior->posterior coupling (downstream)
        w_down = eps;
        w_up = beta*eps;
    end
    W(1,1) = w_m;
    W(1,2) = w_up;
    for ii=2:dim-1
        W(ii, ii) = w_m;
        W(ii,ii-1) = w_down;
        W(ii,ii+1) = w_up;
    end
    W(dim,dim-1) = w_down;
    W(dim,dim) = w_m;

end





end

