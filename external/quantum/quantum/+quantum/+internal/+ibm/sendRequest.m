function output = sendRequest(request, url)
% Internal use only. Sends request to the url endpoint by retrying certain
% failures. This handles all communication with IBM.

% Copyright 2024 The MathWorks, Inc.

maxNumRequests = 3;
for rq = 1:maxNumRequests

    response = send(request, url);

    status = response.StatusCode;
    if status==200
        % OK
        output = response.Body.Data;
        return

    elseif (status==204 && response.Completed)
        % NoContent. This is a special-case for cancelQuantumTask, where it
        % completes but returns NoContent because we send POST request
        % with empty body.
        output = response.Body.Data;
        return

    elseif status==500
        % InternalServerError
        requestDelay = 0.1;
        pause(requestDelay)
        % Retry
        continue

    elseif status==429
        % TooManyRequests. This should return a field indicating how long to
        % wait before making a new request.
        field = getFields(response, "retry-after");
        if ~isempty(field)
            % Safeguard against possible extreme wait time by keeping it
            % below a maximum.
            maxDelay = seconds(hours(1));
            retryDelay = str2double(field.Value);
            requestDelay = min([retryDelay maxDelay]);
        else
            % Delay a few seconds since no retry-after delay was given.
            requestDelay = 5;
        end
        pause(requestDelay)
        % Retry
        continue

    else
        % Retrying the request won't be successful.
        break
    end
end

% Root-cause exception
statusLine = getReasonPhrase(status);
statusCode = double(status);
msg = message("quantum:QuantumDeviceIBM:StatusError", statusCode, statusLine, url);
serverME = MException(msg);

if status==500 || status==429
    % Requests exceed the maximum number of attempts, add exception to
    % indicate the failure is temporary and should be retried.
    baseME = MException(message("quantum:QuantumDeviceIBM:RequestErrorTemporary"));
else
    % The request failed, add exception to indicate request not complete
    baseME = MException(message("quantum:QuantumDeviceIBM:RequestError"));
end

% The user-friendly error message has the specific server failure as its
% cause
throw(addCause(baseME, serverME))

end