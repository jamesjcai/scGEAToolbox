%INPUTSDLG DEMO (Enhanced input dialog box with multiple data types)
%   Demonstrates new "required" field of Formats argument

% Written by: Takeshi Ikuma
% Last Updated: July 17 2013

clear; close all;

Title = 'INPUTSDLG Demo Dialog';

Prompt = {
   'Name','Name','*'
   'Address','Street','*'
   '','Street2',''
   'City','City','*'
   'State','State','*'
   'Zip','Zip','*'
   'Phone','Phone','*'
   'Fax','Fax',''
   '* Required fields','',''
   };
Formats(1,1).required = 'on'; % Name
Formats(1,1).span = [1 3];
Formats(2,1).required = 'on'; % Address
Formats(2,1).span = [1 3];
Formats(3,1).required = 'off'; % Address 2nd line
Formats(3,1).span = [1 3];
Formats(4,1).required = 'on';  % City
Formats(4,2).required = 'on'; % State
Formats(4,3).required = 'on'; % Zip
Formats(5,1).required = 'on'; % phone
Formats(5,2).required = 'off'; %fax
Formats(6,1).type = 'text'; % footnote
Formats(6,1).span = [1 3];

DefAns = struct([]);

Options.AlignControls = 'on';

[Answer,Cancelled] = inputsdlg(Prompt,Title,Formats,DefAns,Options);

