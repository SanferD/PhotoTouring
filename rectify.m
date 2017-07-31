function [] = rectify()
	close all
	im = imread('corridor.jpg');

    [h, w, ~] = size(im);
    w = w*1/4;
    h = h*1/4;
    im = imresize(im, [h w], 'nearest');
       
	f = 4500;
    f = 750;
	px = size(im, 1)/2;
	py = size(im, 2)/2;
	K = [f 0 px; 0 f py; 0 0 1];

	[H, imRect] = RectifyImage(im, K);
    vanish = GetVanishingPoint(imRect);
    
    
 	fTL = [2091 1813 1];
 	fTR = [2671 1815 1];
 	fBR = [2675 2225 1];
 	fBL = [2086 2226 1];
 
 	bTL = [1191 1278 1];
 	bTR = [3507 1286 1];
 	bBR = [3502 2908 1];
 	bBL = [1193 2911 1];

%     fTL = [ 487 483 1 ];
%     fTR = [ 568 481 1 ];
%     fBR = [ 569 539 1 ];
%     fBL = [ 488 539 1 ];
%     
%     bTL = [ 340 402 1 ];
%     bTR = [ 759 403 1 ];
%     bBR = [ 759 711 1 ];
%     bBL = [ 339 712 1 ];
    
	fD = 10;
	bD = GetD(fTL, bTL, vanish, fD, px, py); %2.5

    fwTL = fD*inv(K)*fTL';
	fwTR = fD*inv(K)*fTR';
	fwBR = fD*inv(K)*fBR';
	fwBL = fD*inv(K)*fBL';

	bwTL = bD*inv(K)*bTL';
	bwTR = bD*inv(K)*bTR';
	bwBR = bD*inv(K)*bBR';
	bwBL = bD*inv(K)*bBL';

	fB1 = fwTR - fwTL;
	fB2 = fwBL - fwTL;
	fB3 = fwTL;

	lB1 = bwTL - fwTL;
	lB2 = fwBL - fwTL;
	lB3 = fwTL;

	rB1 = bwTR - fwTR;
	rB2 = fwBR - fwTR;
	rB3 = fwTR;

	uB1 = fwTR - fwTL;
	uB2 = bwTL - fwTL;
	uB3 = fwTL;

	dB1 = fwBL - bwBL;
	dB2 = bwBR - bwBL;
	dB3 = bwBL;

	fH1 = K*[ fB1 fB2 fB3 ];
	lH1 = K*[ lB1 lB2 lB3 ];
	rH1 = K*[ rB1 rB2 rB3 ];
	uH1 = K*[ uB1 uB2 uB3 ];
	dH1 = K*[ dB1 dB2 dB3 ];

	theta = 20/180 * pi;%pi/6;
    theta = -pi/25;
%     theta = 0;
	R = [   cos(theta) 0 sin(theta);
			0 1 0;
			-sin(theta) 0 cos(theta);
		];
	C = [0 0 fD/3];
    C = [0 0 fD/5];
%     C = [0 0 0];
	t = -C';

	fH2 = K*[ R*fB1 R*fB2 (R*fB3 + t) ];
	lH2 = K*[ R*lB1 R*lB2 (R*lB3 + t) ];
	rH2 = K*[ R*rB1 R*rB2 (R*rB3 + t) ];
	uH2 = K*[ R*uB1 R*uB2 (R*uB3 + t) ];
	dH2 = K*[ R*dB1 R*dB2 (R*dB3 + t) ];

	fH = fH2*inv(fH1);
	lH = lH2*inv(lH1);
	rH = rH2*inv(rH1);
	uH = uH2*inv(uH1);
	dH = dH2*inv(dH1);
% return
	[rows, cols, ~] = size(imRect);
	Inside = @(p) p(1)>=1 && p(2)>=1 && p(1)<=rows && p(2)<=cols;
	tour = zeros(size(imRect));

	invfH = inv(fH);
	invlH = inv(lH);
	invrH = inv(rH);
	invuH = inv(uH);
	invdH = inv(dH);

	% myCluster = parcluster('local');
	% myCluster.NumWorkers = 4;  % 'Modified' property now TRUE
	% saveProfile(myCluster);    % 'local' profile now updated,
	% parpool(4);   

	disp('00%')
	for r=1:rows %=3024 y
		for c=1:cols %=4032 x
			coord = [c r 1]';

			% front
			srcCoord = invfH*coord;
			s = round([ srcCoord(1)/srcCoord(3) srcCoord(2)/srcCoord(3) 1 ]);
			a = [ Check(fTL, fTR, s), Check(fTR, fBR, s), Check(fBR, fBL, s), Check(fBL, fTL, s) ];
			if IsValid( NoZeros(a) ) && Inside(s)
			    tour(r, c, :) = imRect(s(2), s(1), :);
			    continue
			end

			% left
			srcCoord = invlH*coord;
			s = round([ srcCoord(1)/srcCoord(3) srcCoord(2)/srcCoord(3) 1 ]);
			a = [ Check(fTL, bTL, s), Check(bTL, bBL, s), Check(bBL, fBL, s), Check(fBL, fTL, s) ];
			if IsValid( NoZeros(a) ) && Inside(s)
			    tour(r, c, :) = imRect(s(2), s(1), :);
			    continue
			end

			% right
  			srcCoord = invrH*coord;
  			s = round([ srcCoord(1)/srcCoord(3) srcCoord(2)/srcCoord(3) 1 ]);
  			a = [ Check(fTR, bTR, s), Check(bTR, bBR, s), Check(bBR, fBR, s), Check(fBR, fTR, s) ];
  			if IsValid( NoZeros(a) ) && Inside(s)
  			    tour(r, c, :) = imRect(s(2), s(1), :);
  			    continue
  			end

			% up
 			srcCoord = invuH*coord;
 			s = round([ srcCoord(1)/srcCoord(3) srcCoord(2)/srcCoord(3) 1 ]);
 			a = [ Check(fTL, fTR, s), Check(fTR, bTR, s), Check(bTR, bTL, s), Check(bTL, fTL, s) ];
 			if IsValid( NoZeros(a) ) && Inside(s)
 			    tour(r, c, :) = imRect(s(2), s(1), :);
 			    continue
 			end

			% down
			srcCoord = invdH*coord;
			s = round([ srcCoord(1)/srcCoord(3) srcCoord(2)/srcCoord(3) 1 ]);
			a = [ Check(fBL, bBL, s), Check(bBL, bBR, s), Check(bBR, fBR, s), Check(fBR, fBL, s) ];
			if IsValid( NoZeros(a) ) && Inside(s)
			    tour(r, c, :) = imRect(s(2), s(1), :);
			    continue
			end

		end
		if mod(r, 25) == 0
			disp(sprintf('\b\b\b\b%0.2d%%', round( (r/rows)*100 )));
		end
	end
	
	out = sprintf('\b\b\b\b100%%');
	disp(out)
	imshow(uint8(tour));

end

function out = NoZeros(X)
	out = [];
	for i=1:numel(X)
		if X(i)~=0
			out = [out X(i)];
		end
	end
end

function valid = IsValid(X)
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

function [H, imRect] = RectifyImage(im, K)
 	pointsTo = [  1251 1671 1;
				  2667 1678 1;
				  3151 2479 1;
				  534 2450 1;
			   ];

%     pointsTo = [    479 488 1;
%                     708 490 1;
%                     787 619 1;
%                     459 615 1;
%                 ];

	pointsFrom = [  0 0 1;
					0 1 1;
					1 1 1;
					1 0 1;
				 ];
			
	homographyGenerator = HomographyGeneratorForImageRectification(K);
	H = homographyGenerator.Generate(pointsFrom, pointsTo);
	imRect = ImageWarping(im, H);
end

function vanish = GetVanishingPoint(im)
% 	l1p1 = [383 798 1];
% 	l1p2 = [781 1033 1];
% 
% 	l2p1 = [1083 2996 1];
% 	l2p2 = [1296 2832 1];
% 
% 	l3p1 = [3637 1204 1];
% 	l3p2 = [3330 1399 1];
% 
% 	l4p1 = [3553 2950 1];
% 	l4p2 = [3338 2776 1];

    l1p1 = [ 196 364 1 ];
    l1p2 = [ 258 389 1 ];
    
    l2p1 = [ 309 748 1 ];
    l2p2 = [ 343 709 1 ];
    
    l3p1 = [ 806 380 1 ];
    l3p2 = [ 755 401 1 ];
    
    l4p1 = [ 763 718 1 ];
    l4p2 = [ 720 679 1 ];

	l1 = getLineFromTwoImagePoints(l1p1, l1p2);
	l2 = getLineFromTwoImagePoints(l2p1, l2p2);
	l3 = getLineFromTwoImagePoints(l3p1, l3p2);
	l4 = getLineFromTwoImagePoints(l4p1, l4p2);

	v = nullSpace([l1; l2; l3; l4]);
	vanish = [v(1)/v(3) v(2)/v(3) v(3)/v(3)];
% 	PlotVanishingLines(im, l1p1, l2p1, l3p1, l4p1, vanish);
end

function line = getLineFromTwoImagePoints( one, two )
	line = cross(one, two);
end

function N = nullSpace(A)
	[~, ~, v] = svd(A);
	N = v(:, end);
end

function bD = GetD(u, v, p, d, px, py)
	A = [ v(1)-px p(1)-px; v(2)-py p(2)-py ];
	b = d*[ u(1)-px; u(2)-py ];
	x = A\b;
	bD = x(1);
end

function [] = PlotVanishingLines(im, one, two, three, four, vanish)
    figure,imshow(im)
    hold on
    drawLine(one, vanish, 'r');
    drawLine(two, vanish, 'r');
    drawLine(three, vanish, 'r');
    drawLine(four, vanish, 'r');
    hold off
end

function [] = drawLine(p1, p2, c)
    line([p1(1) p2(1)], [p1(2) p2(2)], 'Color', c, 'LineWidth', 2)
end

function s = Check(base, tip, src)
    x = tip - base;
    y = src - base;
      
    p = myCross(x, y);
    
    s = sign(p);
end

function p = myCross(x, y)
    p = x(1)*y(2) - x(2)*y(1);
end
