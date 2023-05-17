function opts = parseOpts(def, varargin)
%PARSEOPTS Parse optional inputs for functions into an options structure.
%
%   OPTS = PARSEOPTS(DEF [, STR1, VAL1, ...]) constructs an options
%   structure with field names given by the supplied strings (STR1,
%   STR2 etc.) and corresponding values (VAL1, VAL2, etc.).
%
%   DEF is a cell array of name value pairs defining available options and
%   their respective default values.

% $Id: parseOpts.m,v 1.2 2008-07-25 05:45:51 shaunc Exp $

opts = struct(def{:});

if (length(varargin) == 1)
  error('parseOpts:invalidOptionsList','Invalid options list. Type ''help parseOpts'' for usage information.');
end

for i = 2:2:length(varargin),
  % field names must be strings
  if (~ischar(varargin{i-1})),
    warning('Invalid option name.');
    continue;
  end

  if isfield(opts,varargin{i-1}),
    opts = setfield(opts, varargin{i-1}, varargin{i});
  else
%     error('parseOpts:unknownOption','Unknown option ''%s''.', varargin{i-1});
    warning('Unknown option ''%s''.', varargin{i-1});
  end
  opts = setfield(opts, varargin{i-1}, varargin{i});
end
