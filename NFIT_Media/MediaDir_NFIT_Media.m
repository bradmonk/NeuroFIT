function [Fnames, varargout] = MediaDir_NFIT_Media(varargin)
persistent LookInFolder 
% persistent Files

thisfile='MediaDir_NFIT_Media';
LookInFolder='';

if nargin > 0
	LookInFolder=varargin{1};
end

	% Save current directory
	pw = what;
	cpath = pw.path;
	pd = dir;
	fp = fileparts(cpath);

	% Get full path to this file (AudioMediaDir) or to AudioMediaDir.m
	wfiname = which(thisfile);
	% Save just the path without the file name
	wthisfile = regexprep(wfiname, '(\w*\.m$)', '','once');
	w1 = what(wthisfile);
	
	fulpath = strcat(wthisfile,LookInFolder);
	fulpathslash = strcat(wthisfile,LookInFolder,'*.tif*');
	
	% Simultaniously Save current path and switch path
	% pathNow = cd(w1.path);
	pathNow = cd(fulpath);
	
	% Files = dir;
	% dir

    
	
	
	av_files = dir(fulpathslash);
	
	Nfiles = numel(av_files);

	Fnames = {};
	for nf = 1:Nfiles
		Fnames{nf} = av_files(nf).name;
	end


	% Change path back to currently open folder
	cd(pathNow)

end