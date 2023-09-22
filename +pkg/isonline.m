function yes = isonline(~)

%% Check connection to internet is available
try
    java.net.InetAddress.getByName('google.com');
    yes = true;
catch
    yes = false;
end
end