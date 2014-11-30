function xticksrotated(ax, labels)
	set(ax,'XTick',1:length(labels));
	set(ax,'XTickLabel','');
	% Estimate the location of the labels based on the position 
	% of the xlabel 
	hx = get(ax,'XLabel');  % Handle to xlabel 
	set(hx,'Units','data'); 
	pos = get(hx,'Position'); 
	y = pos(2); 
	% Place the new labels 
	for i = 1:size(labels,1) 
	    t(i) = text(X(i),y,labels(i,:)); 
	end 
	set(t,'Rotation',90,'HorizontalAlignment','right') 
end
