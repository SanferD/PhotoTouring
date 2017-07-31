function [] = tour_photo()
	close all
	im = imread('corridor.jpg');

	% resize the image
    [h, w, ~] = size(im);
    
    % intrinsic params
    f = 4500;
	px = size(im, 1)/2;
	py = size(im, 2)/2;
	K = [f 0 px; 0 f py; 0 0 1];

	Ry = @(angle) [ cos(angle) 0 sin(angle); 0 1 0; -sin(angle) 0 cos(angle); ];
	
	% extrinsic params
	fD = 10;
    thetaStart = -pi/35;
    thetaFinish = -pi/35 - (20/180)*pi;

    CStart = [ 0 0 0 ];
    CFinish = [ 0 0 .2*fD ];

    [Rset, Cset] = InterpolateCoordinate(Ry(thetaStart), CStart, Ry(thetaFinish), CFinish, 10); 

    % compute the tour
	imRect = RectifyImage(im, K);
    vanish = GetVanishingPoint();
	planeGeom = PlaneGeometry(K, fD, vanish);
	[fH1, lH1, rH1, uH1, dH1] = planeGeom.GetFaceHomographies1();
    
	[rows, cols, ~] = size(imRect);
	Inside = @(p) p(1)>=1 && p(2)>=1 && p(1)<=rows && p(2)<=cols;

    for i = 1:numel(Rset)
		% make the tour image
		tour = zeros(size(imRect));
		name = sprintf('output/%d.jpg', i);

		[fH2, lH2, rH2, uH2, dH2] = planeGeom.GetFaceHomographies2(Rset{i}, -Cset(:,i));

		% get the final homography
		fH = fH2*inv(fH1);
		lH = lH2*inv(lH1);
		rH = rH2*inv(rH1);
		uH = uH2*inv(uH1);
		dH = dH2*inv(dH1);

		out = sprintf('%s: 00%%', name);
		disp(out) % display progress
		for r=1:rows 
			for c=1:cols
				coord = [c r 1]';

				s = getSrcCoord(fH, coord);
				if planeGeom.IsInFrontFace(s) && Inside(s)
				    tour(r, c, :) = imRect(s(2), s(1), :);
				    continue
				end

				s = getSrcCoord(lH, coord);
				if planeGeom.IsInLeftFace(s) && Inside(s)
				    tour(r, c, :) = imRect(s(2), s(1), :);
				    continue
				end

				s = getSrcCoord(rH, coord);
				if planeGeom.IsInRightFace(s) && Inside(s)
	  			    tour(r, c, :) = imRect(s(2), s(1), :);
	  			    continue
	  			end

				s = getSrcCoord(uH, coord);
				if planeGeom.IsInUpFace(s) && Inside(s)
	 			    tour(r, c, :) = imRect(s(2), s(1), :);
	 			    continue
	 			end

				s = getSrcCoord(dH, coord);
				if planeGeom.IsInDownFace(s) && Inside(s)
				    tour(r, c, :) = imRect(s(2), s(1), :);
				    continue
				end
			end
			if mod(r, 25) == 0
				disp(sprintf('\b\b\b\b%0.2d%%', round( (r/rows)*100 ))); % update progress
			end
		end
		
		out = sprintf('\b\b\b\b100%%'); % done
		disp(out)

		imwrite(uint8(tour), name);
    end
end

function imRect = RectifyImage(im, K)
 	pointsTo = [  1251 1671 1;
				  2667 1678 1;
				  3151 2479 1;
				  534 2450 1;
			   ];

	pointsFrom = [  0 0 1;
					0 1 1;
					1 1 1;
					1 0 1;
				 ];
			
	homographyGenerator = HomographyGeneratorForImageRectification(K);
	H = homographyGenerator.Generate(pointsFrom, pointsTo);
	imRect = ImageWarping(im, H);
end

function vanish = GetVanishingPoint()
	l1p1 = [383 798 1];
	l1p2 = [781 1033 1];

	l2p1 = [1083 2996 1];
	l2p2 = [1296 2832 1];

	l3p1 = [3637 1204 1];
	l3p2 = [3330 1399 1];

	l4p1 = [3553 2950 1];
	l4p2 = [3338 2776 1];

	l1 = getLineFromTwoImagePoints(l1p1, l1p2);
	l2 = getLineFromTwoImagePoints(l2p1, l2p2);
	l3 = getLineFromTwoImagePoints(l3p1, l3p2);
	l4 = getLineFromTwoImagePoints(l4p1, l4p2);

	v = nullSpace([l1; l2; l3; l4]);
	vanish = [v(1)/v(3) v(2)/v(3) v(3)/v(3)];
end

function line = getLineFromTwoImagePoints( one, two )
	line = cross(one, two);
end

function N = nullSpace(A)
	[~, ~, v] = svd(A);
	N = v(:, end);
end

function s = getSrcCoord(H, c)
	s = H\c;
	s = round( [s(1)/s(3) s(2)/s(3) 1] );
end

function [] = SetupParallel(n)
	myCluster = parcluster('local');
	myCluster.NumWorkers = n;  % 'Modified' property now TRUE
	saveProfile(myCluster);    % 'local' profile now updated,
	parpool(n); 
end
