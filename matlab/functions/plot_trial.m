function plot_trial(trial, fn)
        %plot_trial       Plot the specified solution components as a heatmap for all times specified. Output to gif file
        %              named specified. Use subplot to display multiple solution components at once.
        %
        % Usage:
        %                       plot_trial(trial,fn)
        %
        % Input:
        %                       trial = a trial structure from import_trials
        %                       fn = output file name for plot
        %
        % Examples:
        %                       %plot voltage and A fields
        %                       fn = './test_trial_animation.gif';
        %                       trials = import_trials('Spanky_2013-01-17-1325.mat');
	%			plot_trial(trials(117), fn);

        close all;
	fig = figure('visible', 'off');

        if (nargin < 2)
                throw(MException('Argin:MoreExpected', 'More input arguments expected'));
        end

	nsubplots=3;
	%Plot animation of trial at 60Hz
        frames = 1:length(trial.cursor);
        frames = 1:10;
	framen = 1;
        for i=frames
        	clf(fig);
        	subplot(1,nsubplots,1);
		hold on;
		%Plot cursor data
		plot(trial.cursor(i,1), trial.cursor(i,2),'ro');
		plot(trial.target(1), trial.target(2), 'bo');
		plot(trial.cursorstart(1), trial.cursorstart(2), 'go');
        	xlabel('x');
        	ylabel('y');
		%xlim([-0.5 0.5]); ylim([-0.5 0.5]);
		xlim([-1 1]); ylim([-1 1]);
		title(['Trial. time = ' num2str(trial.times(i))]);
		legend('Cursor', 'Target', 'Start')
		%Plot torque data
		subplot(1,nsubplots, 2);
		plot3(trial.torque(1:i,1), trial.torque(1:i,2), trial.torque(1:i,3))
		xlim([min(trial.torque(:,1)) max(trial.torque(:,1))]); ylim([min(trial.torque(:,2)) max(trial.torque(:,2))]); zlim([min(trial.torque(:,3)) max(trial.torque(:,3))]);
		xlabel('torque axes 1'); ylabel('torque axes 2'); zlabel('torque axes 3');
		title('Torque data')

		%Plot electrode data
		subplot(1,nsubplots, 3)
		cm = hsv(100);
		maxrate4 = max(trial.rates(:,4));
		scatter3(trial.rates(1:i,1), trial.rates(1:i,2), trial.rates(1:i,3), [], cm(floor(100*trial.rates(1:i,4)/maxrate4),:));
		xlim([min(trial.rates(:,1)) max(trial.rates(:,1))]); ylim([min(trial.rates(:,2)) max(trial.rates(:,2))]); zlim([min(trial.rates(:,3)) max(trial.rates(:,3))]);
		xlabel('electrode 1'); ylabel('electrode 2'); zlabel('electrode 3');
		title('Electrode firing rates')

		%Write file to temp png
	        plotmult(gcf, [fn '_tmp.png'], nsubplots, 'png', [6 4]);
	        [X, m] = imread([fn '_tmp.png'], 'png');
 	        f = im2frame(X, m);
         	if i == 1
            		[im, map] = rgb2ind(f.cdata, 256, 'nodither');
            		im(1,1,1,length(frames)) = 0;
            	else
            		im(:,:,1,framen) = rgb2ind(f.cdata, map, 'nodither');
        	end
        	framen = framen + 1;
        end
	%Save to file
	imwrite(im, map, [fn], 'gif', 'DelayTime', 0, 'LoopCount', inf);
end
