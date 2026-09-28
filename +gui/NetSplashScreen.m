classdef NetSplashScreen < handle
    % NetSplashScreen  Frameless splash screen built on Windows .NET forms.
    %
    %   s = gui.NetSplashScreen(imageFile) shows imageFile in a borderless,
    %   always-on-top window centred on the screen, with a thin progress bar
    %   along the bottom edge. It needs no Java runtime, so it works on
    %   MATLAB installations without a bundled JRE; it needs Windows.
    %   Use delete(s) to close it.
    %
    %   It mirrors the part of the gui.SplashScreen (Java/Swing) interface
    %   that gui.sc_splashscreen uses: the ProgressRatio property and the
    %   addText method.
    %
    %   Example:
    %   s = gui.NetSplashScreen("splash.png");
    %   s.addText(30, 50, "SCGEATOOL", FontSize=18, Color=[1 1 1]);
    %   s.ProgressRatio = 0.5;
    %   delete(s)
    %
    %   See also: gui.SplashScreen, gui.sc_splashscreen

    properties
        ProgressRatio (1,1) double {mustBeBetween(ProgressRatio, 0, 1)} = 0  % Fraction of the bar filled
    end

    properties (GetAccess = public, SetAccess = private)
        Width = 0   % Width of the window in pixels
        Height = 0  % Height of the window in pixels
    end

    properties (Access = private)
        Form = []
        Image = []
        Canvas = []
        ProgressTrack = []
        ProgressFill = []
    end

    properties (Constant, Access = private)
        BarMargin = 10    % Gap between the bar and the window's side edges
        BarHeight = 6     % Thickness of the bar
        BarOffset = 5     % Gap between the bar and the window's bottom edge
    end

    methods

        function obj = NetSplashScreen(imageFile)
            % Construct the splash screen and put it on screen.
            arguments
                imageFile (1,1) string {mustBeFile}
            end
            NET.addAssembly("System.Windows.Forms");
            NET.addAssembly("System.Drawing");

            obj.Image = System.Drawing.Image.FromFile(char(imageFile));
            obj.Width = obj.Image.Width;
            obj.Height = obj.Image.Height;

            form = System.Windows.Forms.Form();
            form.FormBorderStyle = System.Windows.Forms.FormBorderStyle.None;
            form.StartPosition = System.Windows.Forms.FormStartPosition.CenterScreen;
            form.TopMost = true;
            form.ShowInTaskbar = false;
            form.ClientSize = System.Drawing.Size(obj.Width, obj.Height);
            % Text is drawn into this copy, as the Java splash does. Label
            % controls cannot be used: a transparent Label shows only the
            % form's background through it, so overlapping labels clip
            % each other.
            obj.Canvas = System.Drawing.Bitmap(obj.Image);
            form.BackgroundImage = obj.Canvas;
            form.BackgroundImageLayout = System.Windows.Forms.ImageLayout.None;
            obj.Form = form;

            % Two plain panels rather than a ProgressBar control, to match
            % the thin grey bar of the Java splash.
            barWidth = obj.Width - 2*obj.BarMargin;
            barTop = obj.Height - obj.BarOffset - obj.BarHeight;
            obj.ProgressTrack = System.Windows.Forms.Panel();
            obj.ProgressTrack.SetBounds(obj.BarMargin, barTop, barWidth, obj.BarHeight);
            obj.ProgressTrack.BackColor = System.Drawing.Color.FromArgb(int32(64), int32(64), int32(64));
            obj.ProgressFill = System.Windows.Forms.Panel();
            obj.ProgressFill.SetBounds(0, 0, 0, obj.BarHeight);
            obj.ProgressFill.BackColor = System.Drawing.Color.FromArgb(int32(190), int32(190), int32(190));
            obj.ProgressTrack.Controls.Add(obj.ProgressFill);
            form.Controls.Add(obj.ProgressTrack);

            form.Show();
            obj.refresh();
        end

        function addText(obj, x, y, text, options)
            % addText  Write text on the splash screen.
            %
            %   S.addText(X,Y,STR) writes STR with its baseline starting at
            %   pixel [X,Y], the same convention gui.SplashScreen uses.
            %   Name-value options: FontSize (points, at the nominal 96
            %   pixels per inch of the image), FontName, Color (an RGB
            %   triplet or any color name validatecolor accepts), and
            %   Shadow (true by default).
            arguments
                obj
                x (1,1) double
                y (1,1) double
                text (1,1) string
                options.FontSize (1,1) double {mustBePositive} = 12
                options.FontName (1,1) string = "Arial"
                options.Color = [1 0 0]
                options.Shadow (1,1) logical = true
            end
            rgb = round(255*validatecolor(options.Color));
            % Size the font in image pixels so the layout does not change
            % with the display's scaling.
            pixelsPerPoint = 96/72;
            fontPixels = options.FontSize*pixelsPerPoint;
            style = System.Drawing.FontStyle.Regular;
            font = System.Drawing.Font(char(options.FontName), single(fontPixels), ...
                style, System.Drawing.GraphicsUnit.Pixel);
            family = font.FontFamily;
            ascent = fontPixels*double(family.GetCellAscent(style))/double(family.GetEmHeight(style));

            gfx = System.Drawing.Graphics.FromImage(obj.Canvas);
            gfx.TextRenderingHint = System.Drawing.Text.TextRenderingHint.AntiAlias;
            % GenericTypographic drops the padding DrawString adds, so the
            % top-left of the layout box is exactly one ascent above the
            % baseline.
            format = System.Drawing.StringFormat.GenericTypographic;
            top = y - ascent;
            if options.Shadow
                shadowAlpha = int32(110);
                obj.drawString(gfx, text, font, System.Drawing.Color.FromArgb(shadowAlpha, ...
                    int32(0), int32(0), int32(0)), x + 1, top + 1, format);
            end
            obj.drawString(gfx, text, font, System.Drawing.Color.FromArgb(int32(rgb(1)), ...
                int32(rgb(2)), int32(rgb(3))), x, top, format);
            gfx.Dispose();
            font.Dispose();
            obj.Form.Invalidate();
            obj.refresh();
        end

        function set.ProgressRatio(obj, value)
            obj.ProgressRatio = value;
            obj.updateProgressBar();
        end

        function delete(obj)
            % delete  Close the window and release the image file.
            if ~isempty(obj.Form)
                obj.Form.Close();
                obj.Form.Dispose();
            end
            if ~isempty(obj.Canvas)
                obj.Canvas.Dispose();
            end
            if ~isempty(obj.Image)
                obj.Image.Dispose();
            end
        end

    end

    methods (Access = private)

        function updateProgressBar(obj)
            if isempty(obj.ProgressFill)
                return
            end
            obj.ProgressFill.Width = round(obj.ProgressRatio*obj.ProgressTrack.Width);
            obj.refresh();
        end

        function drawString(~, gfx, text, font, color, x, y, format)
            brush = System.Drawing.SolidBrush(color);
            gfx.DrawString(char(text), font, brush, single(x), single(y), format);
            brush.Dispose();
        end

        function refresh(obj)
            % Nothing runs a message loop for this window while the app
            % starts up, so repaint it explicitly after every change.
            obj.Form.Refresh();
            System.Windows.Forms.Application.DoEvents();
        end

    end

end
