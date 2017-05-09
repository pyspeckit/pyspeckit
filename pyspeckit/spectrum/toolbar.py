from __future__ import print_function
from matplotlib.backend_bases import NavigationToolbar2

class NavigationToolbar3(NavigationToolbar2):

    def press_pan(self, event):
        print("pan pressed")
        super(NavigationToolbar3,self).press_pan(self,event)

class MyNavToolbar(NavigationToolbar2):
    """wx/mpl NavToolbar hack with an additional tools user interaction.
    This class is necessary because simply adding a new togglable tool to the
    toolbar won't (1) radio-toggle between the new tool and the pan/zoom tools.
    (2) disable the pan/zoom tool modes in the associated subplot(s).
    """
    def __init__(self, canvas):
        super(MyNavToolbar, self).__init__(canvas)
        self.pan_tool  = self.FindById(self._NTB2_PAN)
        self.zoom_tool = self.FindById(self._NTB2_ZOOM)
        self.Bind(wx.EVT_TOOL, self.on_toggle_pan_zoom, self.zoom_tool)
        self.Bind(wx.EVT_TOOL, self.on_toggle_pan_zoom, self.pan_tool)

        self.user_tools = {}   # user_tools['tool_mode'] : wx.ToolBarToolBase

        self.InsertSeparator(5)
        self.add_user_tool('lasso', 6, icons.lasso_tool.ConvertToBitmap(), True, 'Lasso')
        self.add_user_tool('gate', 7, icons.gate_tool.ConvertToBitmap(), True, 'Gate')

    def add_user_tool(self, mode, pos, bmp, istoggle=True, shortHelp=''):
        """Adds a new user-defined tool to the toolbar.
        mode -- the value that MyNavToolbar.get_mode() will return if this tool 
                is toggled on
        pos -- the position in the toolbar to add the icon
        bmp -- a wx.Bitmap of the icon to use in the toolbar
        isToggle -- whether or not the new tool toggles on/off with the other 
                    togglable tools
        shortHelp -- the tooltip shown to the user for the new tool
        """
        tool_id = wx.NewId()
        self.user_tools[mode] = self.InsertSimpleTool(pos, tool_id, bmp,
                            isToggle=istoggle, shortHelpString=shortHelp)
        self.Bind(wx.EVT_TOOL, self.on_toggle_user_tool, self.user_tools[mode])

    def get_mode(self):
        """Use this rather than navtoolbar.mode
        """
        for mode, tool in self.user_tools.items():
            if tool.IsToggled():
                return mode
        return self.mode

    def untoggle_mpl_tools(self):
        """Hack city: Since I can't figure out how to change the way the 
        associated subplot(s) handles mouse events: I generate events to turn
        off whichever tool mode is enabled (if any). 
        This function needs to be called whenever any user-defined tool 
        (eg: lasso) is clicked.
        """
        if self.pan_tool.IsToggled():
            wx.PostEvent(
                self.GetEventHandler(), 
                wx.CommandEvent(wx.EVT_TOOL.typeId, self._NTB2_PAN)
            )
            self.ToggleTool(self._NTB2_PAN, False)
        elif self.zoom_tool.IsToggled():
            wx.PostEvent(
                self.GetEventHandler(),
                wx.CommandEvent(wx.EVT_TOOL.typeId, self._NTB2_ZOOM)
            )
            self.ToggleTool(self._NTB2_ZOOM, False)

    def on_toggle_user_tool(self, evt):
        """user tool click handler.
        """
        if evt.Checked():
            self.untoggle_mpl_tools()
            #untoggle other user tools
            for tool in self.user_tools.values():
                if tool.Id != evt.Id:
                    self.ToggleTool(tool.Id, False)

    def on_toggle_pan_zoom(self, evt):
        """Called when pan or zoom is toggled. 
        We need to manually untoggle user-defined tools.
        """
        if evt.Checked():
            for tool in self.user_tools.values():
                self.ToggleTool(tool.Id, False)
        # Make sure the regular pan/zoom handlers get the event
        evt.Skip()

    def reset_history(self):
        """More hacky junk to clear/reset the toolbar history.
        """
        self._views.clear()
        self._positions.clear()
        self.push_current()

