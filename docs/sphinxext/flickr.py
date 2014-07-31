from docutils import nodes
from docutils.parsers.rst import directives

CODE = """\

<object width="%(width)i" height="%(height)i"> 
<param name="flashvars" value="offsite=true&lang=en-us&page_show_url=%2Fphotos%2F75341362%40N04%2Fsets%2F72157629059614751%2Fshow%2F&page_show_back_url=%2Fphotos%2F75341362%40N04%2Fsets%2F72157629059614751%2F&set_id=72157629059614751&jump_to="> </param>
<param name="movie" value="http://www.flickr.com/apps/slideshow/show.swf?v=%(flickid)s"> </param>
<param name="allowFullScreen" value="true"></param>
<embed type="application/x-shockwave-flash" src="http://www.flickr.com/apps/slideshow/show.swf?v=%(flickid)s" allowFullScreen="true" flashvars="offsite=true&lang=en-us&page_show_url=%2Fphotos%2F75341362%40N04%2Fsets%2F72157629059614751%2Fshow%2F&page_show_back_url=%2Fphotos%2F75341362%40N04%2Fsets%2F72157629059614751%2F&set_id=72157629059614751&jump_to=" width="%(width)i" height="%(height)i"></embed>
</object>
"""

PARAM = """\n    <param name="%s" value="%s"></param>"""

def flickr(name, args, options, content, lineno,
            contentOffset, blockText, state, stateMachine):
    """ Restructured text extension for inserting flickr embedded slideshows """
    if len(content) == 0:
        return
    string_vars = {
        'flickid': content[0],
        'width': 400,
        'height': 300,
        'extra': ''
        }
    extra_args = content[1:] # Because content[0] is ID
    extra_args = [ea.strip().split("=") for ea in extra_args] # key=value
    extra_args = [ea for ea in extra_args if len(ea) == 2] # drop bad lines
    extra_args = dict(extra_args)
    if 'width' in extra_args:
        string_vars['width'] = extra_args.pop('width')
    if 'height' in extra_args:
        string_vars['height'] = extra_args.pop('height')
    if extra_args:
        params = [PARAM % (key, extra_args[key]) for key in extra_args]
        string_vars['extra'] = "".join(params)
    return [nodes.raw('', CODE % (string_vars), format='html')]
flickr.content = True
directives.register_directive('flickr', flickr)

from sphinx.util.compat import Directive

class FlickrDirective(Directive):

    def run(self, *args):
        return [flickr(*args)]

def setup(app):
    app.add_config_value('flickr',False,False)
    app.add_config_value('flickrID', None, 'html')
    app.add_directive('flickr',FlickrDirective)
