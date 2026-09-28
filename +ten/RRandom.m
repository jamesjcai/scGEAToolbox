classdef RRandom < handle
    %RRANDOM  R's default random number generator.
    %   STREAM = TEN.RRANDOM(SEED) is the equivalent of set.seed(SEED) in R with
    %   R's defaults (Mersenne-Twister, Inversion, Rejection). Its methods
    %   return the same values as R's runif, rnorm and sample, so the
    %   reference mode of TEN.SCTENIFOLDNET and TEN.SCTENIFOLDKNK draws the
    %   same cells and the same initial CP factors as the R packages.
    %
    %   MATLAB's own 'twister' cannot stand in for it: it is the same
    %   generator, but seeds it differently and builds each double from two
    %   32-bit draws, where R uses one.
    %
    %   Methods:
    %     u = stream.unifRand(n)          % n draws, as runif(n)
    %     z = stream.rnorm(n)             % n draws, as rnorm(n)
    %     i = stream.sample(n, size)      % as sample(n, size, replace = TRUE)
    %
    %   Example:
    %     stream = ten.RRandom(1);
    %     stream.sample(1000, 8)           % 836 679 129 930 509 471 299 270
    %
    %   Port of scTenifoldpy's RRandom (scTenifold/core/_rng.py), itself
    %   checked against R 4.5.
    %
    % see also: TEN.SCTENIFOLDNET, TEN.SCTENIFOLDKNK, TEN.I_CPALS

    properties (Access = private)
        State (624, 1) uint32
        Position (1, 1) double
    end

    properties (Constant, Access = private)
        N = 624
        M = 397
        MatrixA = uint32(0x9908b0df)
        UpperMask = uint32(0x80000000)
        LowerMask = uint32(0x7fffffff)
        TwoPow32Inv = 2.3283064365386963e-10    % 1 / 2^32, as MT_genrand
        HalfUlp = 0.5 * 2.328306437080797e-10   % 0.5 / (2^32 - 1), as fixup
        Big = 134217728                         % 2^27, as rnorm's inversion
    end

    methods
        function obj = RRandom(seed)
            arguments
                seed (1, 1) double {mustBeInteger} = 1
            end
            % set.seed(): scramble the seed with 50 LCG steps, then one step
            % per element of .Random.seed. Element 1 is the position in the
            % state, overwritten with 624 so that the first draw regenerates.
            s = mod(seed, 2^32);
            for j = 1:50
                s = mod(69069*s + 1, 2^32);
            end
            key = zeros(obj.N + 1, 1);
            for j = 1:(obj.N + 1)
                s = mod(69069*s + 1, 2^32);
                key(j) = s;
            end
            obj.State = uint32(key(2:end));
            obj.Position = obj.N;
        end

        function u = unifRand(obj, n)
            %UNIFRAND  N uniform draws in (0, 1), as R's unif_rand.
            u = double(obj.nextInt32(n))*obj.TwoPow32Inv;
            u(u <= 0) = obj.HalfUlp;
            u((1 - u) <= 0) = 1 - obj.HalfUlp;
        end

        function z = rnorm(obj, n)
            %RNORM  N standard normal draws, as R's rnorm(n).
            % Inversion with two uniforms per draw, as R's norm_rand
            u = obj.unifRand(2*n);
            u = floor(obj.Big*u(1:2:end)) + u(2:2:end);
            z = ten.RRandom.qnorm(u/obj.Big);
        end

        function idx = sample(obj, n, sampleSize)
            %SAMPLE  SAMPLESIZE draws from 1:N with replacement.
            %   Same as sample(N, SAMPLESIZE, replace = TRUE) in R.
            idx = zeros(sampleSize, 1);
            for k = 1:sampleSize
                idx(k) = obj.unifIndex(n) + 1;
            end
        end
    end

    methods (Access = private)
        function v = unifIndex(obj, n)
            % R_unif_index(): rejection sampling from the integers below the
            % next power of two, built 16 bits per uniform (rbits)
            if n <= 0
                v = 0;
                return
            end
            bits = ceil(log2(n));
            nChunks = floor(bits/16) + 1;
            while true
                u = obj.unifRand(nChunks);
                v = 0;
                for k = 1:nChunks
                    v = 65536*v + floor(u(k)*65536);
                end
                v = mod(v, 2^bits);
                if v < n
                    return
                end
            end
        end

        function y = nextInt32(obj, n)
            y = zeros(n, 1, "uint32");
            for k = 1:n
                if obj.Position >= obj.N
                    obj.twist();
                end
                obj.Position = obj.Position + 1;
                y(k) = obj.State(obj.Position);
            end
            % Tempering
            y = bitxor(y, bitshift(y, -11));
            y = bitxor(y, bitand(bitshift(y, 7), uint32(0x9d2c5680)));
            y = bitxor(y, bitand(bitshift(y, 15), uint32(0xefc60000)));
            y = bitxor(y, bitshift(y, -18));
        end

        function twist(obj)
            % MT19937 state regeneration. Each block reads only entries that
            % the sequential loop would already have updated (or not yet
            % touched), so it can be vectorized.
            mt = obj.State;
            blocks = {1:227, 228:454, 455:623};
            for b = 1:numel(blocks)
                kk = blocks{b};
                mt(kk) = obj.mix(mt(kk), mt(kk + 1), mt(mod(kk + obj.M - 1, obj.N) + 1));
            end
            mt(obj.N) = obj.mix(mt(obj.N), mt(1), mt(obj.M));
            obj.State = mt;
            obj.Position = 0;
        end

        function out = mix(obj, current, next, far)
            y = bitor(bitand(current, obj.UpperMask), bitand(next, obj.LowerMask));
            out = bitxor(bitxor(far, bitshift(y, -1)), obj.MatrixA*bitand(y, uint32(1)));
        end
    end

    methods (Static)
        function x = qnorm(p)
            %QNORM  Standard normal quantile, as R's qnorm (Wichura's AS 241).
            q = p - 0.5;
            x = zeros(size(p));

            central = abs(q) <= 0.425;
            qc = q(central);
            r = 0.180625 - qc.*qc;
            x(central) = qc.*(((((((r*2509.0809287301226727 + ...
                33430.575583588128105).*r + 67265.770927008700853).*r + ...
                45921.953931549871457).*r + 13731.693765509461125).*r + ...
                1971.5909503065514427).*r + 133.14166789178437745).*r + ...
                3.387132872796366608) ...
                ./(((((((r*5226.495278852545925 + ...
                28729.085735721942674).*r + 39307.89580009271061).*r + ...
                21213.794301586595867).*r + 5394.1960214247511077).*r + ...
                687.1870074920579083).*r + 42.313330701600911252).*r + 1.0);

            tail = ~central;
            qt = q(tail);
            pt = p(tail);
            pSmall = pt;
            pSmall(qt > 0) = 1 - pt(qt > 0);
            r = sqrt(-log(pSmall));
            val = zeros(size(r));
            near = r <= 5.0;
            rn = r(near) - 1.6;
            val(near) = (((((((rn*7.7454501427834140764e-4 + ...
                0.0227238449892691845833).*rn + 0.24178072517745061177).*rn + ...
                1.27045825245236838258).*rn + 3.64784832476320460504).*rn + ...
                5.7694972214606914055).*rn + 4.6303378461565452959).*rn + ...
                1.42343711074968357734) ...
                ./(((((((rn*1.05075007164441684324e-9 + ...
                5.475938084995344946e-4).*rn + 0.0151986665636164571966).*rn + ...
                0.14810397642748007459).*rn + 0.68976733498510000455).*rn + ...
                1.6763848301838038494).*rn + 2.05319162663775882187).*rn + 1.0);
            rf = r(~near) - 5.0;
            val(~near) = (((((((rf*2.01033439929228813265e-7 + ...
                2.71155556874348757815e-5).*rf + 0.0012426609473880784386).*rf + ...
                0.026532189526576123093).*rf + 0.29656057182850489123).*rf + ...
                1.7848265399172913358).*rf + 5.4637849111641143699).*rf + ...
                6.6579046435011037772) ...
                ./(((((((rf*2.04426310338993978564e-15 + ...
                1.4215117583164458887e-7).*rf + 1.8463183175100546818e-5).*rf + ...
                7.868691311456132591e-4).*rf + 0.0148753612908506148525).*rf + ...
                0.13692988092273580531).*rf + 0.59983220655588793769).*rf + 1.0);
            val(qt < 0) = -val(qt < 0);
            x(tail) = val;
        end
    end
end
