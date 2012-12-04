function p = diff_vec(u,v)
% p = diff_vec(u,v)
% Returns a difference ratio of v, w.r.t. u.
p = norm(u - v) / norm(u);
