function [T] = virtualknockout(obj, gid)
if nargin < 2, gid = 1; end
assert(floor(gid) == gid);
assert(gid >= 1 & gid <= obj.NumGenes);
import ten.*
[T] = i_knk(obj.A, gid, obj.g);
