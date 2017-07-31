classdef HomographyGeneratorForImageRectification < handle

	properties (Access = private)
		K
	end

	methods
		function obj = HomographyGeneratorForImageRectification(K)
			obj.K = K;
		end

		function HForRectification = Generate(obj, pointsFrom, pointsTo)
			HToGroundPlane = ComputeHomography(pointsTo, pointsFrom);
			HForRectification = obj.ComputeHomographyForRectification(HToGroundPlane);
		end
	end

	methods (Access = private)

		function H_new = ComputeHomographyForRectification(obj, H)
			RCam = obj.GetRotationFromHomography(H);
			RCamToRRect = obj.GetRotationToGroundPlane(RCam);
			H_new = obj.GetHomographyFromTwoRotations(RCam, RCamToRRect);
		end

		function R = GetRotationFromHomography(obj, H)
			HomographyManipulator = HomographyFromImageToGroundPlaneManipulator(obj.K, H);
			R = HomographyManipulator.GetRotation();
		end

		function RNew = GetRotationToGroundPlane(obj, R)
			cameraYAxis = [0 0 -1];
			cameraXAxis = obj.GetCameraXAxis(R(1,:), cameraYAxis);
			cameraZAxis = obj.GetCameraZAxisFromXYAxis(cameraXAxis, cameraYAxis);			
			RNew = [cameraXAxis; cameraYAxis; cameraZAxis];
		end

		function newX = GetCameraXAxis(obj, oldX, newY);
			oldXonNewY = obj.ProjectVectorOntoAxis(oldX, newY);
			newX = oldX - oldXonNewY;
		end	

		function projection = ProjectVectorOntoAxis(obj, vec, axis)
			projection = dot(vec, axis)*axis/norm(axis);
		end

		function r3 = GetCameraZAxisFromXYAxis(obj, r1, r2)
			r3 = cross(r1,r2);
		end

		function HNew = GetHomographyFromTwoRotations(obj, R, RNew)
			HNew = obj.K * RNew * inv(R) * inv(obj.K);
		end

	end
end
