% Rank Biased Overlap (RBO)

https: /  / towardsdatascience.com / rbo - v - s - kendall - tau - to - compare - ranked - lists - of - items - 8776c5182899

https: /  / stackoverflow.com / questions / 13574406 / how - to - compare - ranked - lists

https: /  / ai.plainenglish.io / comparing - top - k - rankings - statistically - 9adfc9cfc98b

https: /  / github.com / dlukes / rbo / blob / master / rbo.py

https: /  / github.com / changyaochen / rbo

In [1]:import rbo

In [2]:S = [1, 2, 3]

In [3]:T = [1, 3, 2]

In [4]:rbo.RankingSimilarity(S, T).rbo()
Out[4]:0.8333333333333334


% s=py.list(1,2,3);
% t=py.list(1,3,2);
% rob = py.importlib.import_module('rbo');
% rob.RankingSimilarity(s,t).rob()


pyrun("import rbo")
S = pyrun("S = [1, 2, 3]", "S");
T = pyrun("T = [1, 3, 2]", "T");
%pyrun("S = s", s=[1, 2, 3]);
%pyrun("T = t", t=[1, 3, 2]);
r = pyrun("r=rbo.RankingSimilarity(S, T).rbo()", "r");
