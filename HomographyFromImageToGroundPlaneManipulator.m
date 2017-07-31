classdef HomographyFromImageToGroundPlaneManipulator < handle

	properties (Access = private)
		invK;
		row1Norm;
		H;
	end

	methods 
		function obj = HomographyFromImageToGroundPlaneManipulator(K, H)
			obj.invK = inv(K);
			obj.H = H;
			obj.row1Norm = norm(obj.invK*H(:,1));
		end

		function R = GetRotation(obj)
			r1 = obj.GetFirstRotationRow();
			r2 = obj.GetSecondRotationRow();
			r3 = obj.GetThirdRotationRowFromThePreviousTwo(r1, r2);
			R = [r1 r2 r3];
		end
	end

	methods (Access = private)
		
		function r1 = GetFirstRotationRow(obj)
			r1 = obj.GetExtrinsicParamFromHomography(obj.H(:, 1));
		end

		function r2 = GetSecondRotationRow(obj)
			r2 = obj.GetExtrinsicParamFromHomography(obj.H(:, 2));
		end

		function r3 = GetThirdRotationRowFromThePreviousTwo(obj, r1, r2)
			r3 = cross(r1, r2);
		end

		function t = GetTranslation(obj)
			t = obj.GetExtrinsicParamFromHomography(obj.H(:, 3));
		end

		% Since H = K*[r1 r2 t], [r1 r2 t] = [inv(K)*h1 inv(K)*h2 inv(K)*t];
		% We use the scale factor norm(r1) so that each rotation column is ~1.
		function extrinsicCol = GetExtrinsicParamFromHomography(obj, homographyCol)
			extrinsicCol = obj.invK*homographyCol;
			extrinsicCol = obj.DivideByNormOfRow1ToNormalize(extrinsicCol);
		end

		function n = DivideByNormOfRow1ToNormalize(obj, x)
			n = x/obj.row1Norm;
		end

	end

end