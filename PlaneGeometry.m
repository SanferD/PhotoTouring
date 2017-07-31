classdef PlaneGeometry < handle

	properties (Access = private)
		fTL; fTR; fBR; fBL;
		bTL; bTR; bBR; bBL;

		fB1; fB2; fB3;
		lB1; lB2; lB3;
		rB1; rB2; rB3;
		uB1; uB2; uB3;
		dB1; dB2; dB3;

		K;
	end

	methods (Access = public)
		function obj = PlaneGeometry(K, fd, vanish)
		    % f = front
		 	obj.fTL = [2091 1813 1];
		 	obj.fTR = [2671 1815 1];
		 	obj.fBR = [2675 2225 1];
		 	obj.fBL = [2086 2226 1];
		 
		 	obj.bTL = [1191 1278 1];
		 	obj.bTR = [3507 1286 1];
		 	obj.bBR = [3502 2908 1];
		 	obj.bBL = [1193 2911 1];

		    % obj.fTL = [ 487 483 1 ]; % TL = Top Left
		    % obj.fTR = [ 568 481 1 ]; % TR = Top Right
		    % obj.fBR = [ 569 539 1 ]; % BR = Back Right
		    % obj.fBL = [ 488 539 1 ]; % BL = Back Left

		    % % b = back
		    % obj.bTL = [ 340 402 1 ];
		    % obj.bTR = [ 759 403 1 ];
		    % obj.bBR = [ 759 711 1 ];
		    % obj.bBL = [ 339 712 1 ];

			bd = obj.getD(obj.fTL, obj.bTL, vanish, fd, K(1,3), K(2,3)); %2.5
		    
		    % world coordinates
		    iK = inv(K);
		    computeFrontWorldCoordinates = @(v) obj.computeWorldCoordinates(iK, fd, v');
		    computeBackWorldCoordinates = @(v) obj.computeWorldCoordinates(iK, bd, v');

		    fwTL = computeFrontWorldCoordinates(obj.fTL);
			fwTR = computeFrontWorldCoordinates(obj.fTR);
			fwBR = computeFrontWorldCoordinates(obj.fBR);
			fwBL = computeFrontWorldCoordinates(obj.fBL);

			bwTL = computeBackWorldCoordinates(obj.bTL);
			bwTR = computeBackWorldCoordinates(obj.bTR);
			bwBR = computeBackWorldCoordinates(obj.bBR);
			bwBL = computeBackWorldCoordinates(obj.bBL);

			% basis vectors
			[obj.fB1, obj.fB2, obj.fB3] = obj.computeBasisVectors( fwTL, fwTR, fwBL );
			[obj.lB1, obj.lB2, obj.lB3] = obj.computeBasisVectors( fwTL, bwTL, fwBL );
			[obj.rB1, obj.rB2, obj.rB3] = obj.computeBasisVectors( fwTR, bwTR, fwBR );
			[obj.uB1, obj.uB2, obj.uB3] = obj.computeBasisVectors( fwTL, fwTR, bwTL );
			[obj.dB1, obj.dB2, obj.dB3] = obj.computeBasisVectors( bwBL, fwBL, bwBR );
		
			obj.K = K;
		end

		function [f, l, r, u, d] = GetFaceHomographies1(obj)
			f = obj.getFrontHomography1();
			l = obj.getLeftHomography1();
			r = obj.getRightHomography1();
			u = obj.getUpHomography1();
			d = obj.getDownHomography1();
		end

		function [f, l, r, u, d] = GetFaceHomographies2(obj, R, T)
			f = obj.getFrontHomography2(R, T);
			l = obj.getLeftHomography2(R, T);
			r = obj.getRightHomography2(R, T);
			u = obj.getUpHomography2(R, T);
			d = obj.getDownHomography2(R, T);
		end

		function b = IsInFrontFace(obj, s)
			one = obj.getSign(obj.fTL, obj.fTR, s);
			two = obj.getSign(obj.fTR, obj.fBR, s);
			three = obj.getSign(obj.fBR, obj.fBL, s);
			four = obj.getSign(obj.fBL, obj.fTL, s);
			b = obj.isValid( obj.removeZeros([ one, two, three, four ]) );
		end

		function b = IsInLeftFace(obj, s)
			one = obj.getSign(obj.fTL, obj.bTL, s);
			two = obj.getSign(obj.bTL, obj.bBL, s);
			three = obj.getSign(obj.bBL, obj.fBL, s);
			four = obj.getSign(obj.fBL, obj.fTL, s);
			b = obj.isValid( obj.removeZeros([ one, two, three, four ]) );
		end

		function b = IsInRightFace(obj, s)
			one = obj.getSign(obj.fTR, obj.bTR, s);
			two = obj.getSign(obj.bTR, obj.bBR, s);
			three = obj.getSign(obj.bBR, obj.fBR, s);
			four = obj.getSign(obj.fBR, obj.fTR, s);
			b = obj.isValid( obj.removeZeros([ one, two, three, four ]) );
		end

		function b = IsInUpFace(obj, s)
			one = obj.getSign(obj.fTL, obj.fTR, s);
			two = obj.getSign(obj.fTR, obj.bTR, s);
			three = obj.getSign(obj.bTR, obj.bTL, s);
			four = obj.getSign(obj.bTL, obj.fTL, s);
			b = obj.isValid( obj.removeZeros([ one, two, three, four ]) );
		end

		function b = IsInDownFace(obj, s)
			one = obj.getSign(obj.fBL, obj.bBL, s);
			two = obj.getSign(obj.bBL, obj.bBR, s);
			three = obj.getSign(obj.bBR, obj.fBR, s);
			four = obj.getSign(obj.fBR, obj.fBL, s);
			b = obj.isValid( obj.removeZeros([ one, two, three, four ]) );
		end
	end

	methods (Access = private)
		function bD = getD(obj, u, v, p, d, px, py)
			A = [ v(1)-px p(1)-px; v(2)-py p(2)-py ];
			b = d*[ u(1)-px; u(2)-py ];
			x = A\b;
			bD = x(1);
		end

		function w = computeWorldCoordinates(obj, K, d, v)
			w = d*K*v;
		end

		function [B1, B2, B3] = computeBasisVectors(obj, base, tip1, tip2)
			B1 = tip1 - base;
			B2 = tip2 - base;
			B3 = base;
		end

		function H = getFrontHomography1(obj)
			H = obj.K*[ obj.fB1 obj.fB2 obj.fB3 ];
		end

		function H = getLeftHomography1(obj)
			H = obj.K*[ obj.lB1 obj.lB2 obj.lB3 ];
		end

		function H = getRightHomography1(obj)
			H = obj.K*[ obj.rB1 obj.rB2 obj.rB3 ];
		end

		function H = getUpHomography1(obj)
			H = obj.K*[ obj.uB1 obj.uB2 obj.uB3 ];
		end

		function H = getDownHomography1(obj)
			H = obj.K*[ obj.dB1 obj.dB2 obj.dB3 ];
		end
		
		function H = getFrontHomography2(obj, R, t)
			H = obj.K*[ R*obj.fB1 R*obj.fB2 (R*obj.fB3 + t) ];
		end

		function H = getLeftHomography2(obj, R, t)
			H = obj.K*[ R*obj.lB1 R*obj.lB2 (R*obj.lB3 + t) ];
		end

		function H = getRightHomography2(obj, R, t)
			H = obj.K*[ R*obj.rB1 R*obj.rB2 (R*obj.rB3 + t) ];
		end

		function H = getUpHomography2(obj, R, t)
			H = obj.K*[ R*obj.uB1 R*obj.uB2 (R*obj.uB3 + t) ];
		end

		function H = getDownHomography2(obj, R, t)
			H = obj.K*[ R*obj.dB1 R*obj.dB2 (R*obj.dB3 + t) ];
		end
	
		function s = getSign(obj, base, tip, src)
		    x = tip - base;
		    y = src - base;
		    p = obj.myCrossProduct(x, y);
		    s = sign(p);
		end

		function p = myCrossProduct(obj, x, y)
    		p = x(1)*y(2) - x(2)*y(1);
		end

		function out = removeZeros(obj, X)
			out = [];
			for i=1:numel(X)
				if X(i)~=0
					out = [out X(i)];
				end
			end
		end

		function valid = isValid(obj, X)
			if isempty(X)
				valid = false;
				return
			end

			s = sign(X(1));
			for i=2:numel(X)
				if s~=sign(X(i))
					valid = false;
					return
				end
			end

			valid = true;
		end
	end

end