class Font:
    def __init__(self, face='helvetica', size=10, bold=False, italic=False):
        if face not in ('times', 'helvetica', 'courier'):
            face = 'helvetica'
        self.face = face
        self.size = int(size)
        self.bold = bold
        self.italic = italic

    def __str__(self):
        bold = 'bold' if self.bold else 'normal'
        italic = 'italic' if self.italic else 'roman'
        return f'"{self.face} {self.size} {bold} {italic}"'

    def __eq__(self, other):

        eq = (isinstance(self, type(other))) and \
            (self.face == other.face) and \
            (self.size == other.size) and \
            (self.bold == other.bold) and \
            (self.italic == other.italic)
        return eq


class Attributes(dict):
    def __init__(self, color='green', dashlist='8 3', font=None, select=True,
                 highlite=True, dash=False, fixed=False, edit=True, move=True,
                 rotate=True, delete=True, include=True, source=True, text=None,
                 width=1):

        # set the attributes for the region
        self['color'] = color
        self['dashlist'] = dashlist
        if not isinstance(font, Font):
            if isinstance(font, (list, tuple)):
                font = Font(*font)
            elif isinstance(font, dict):
                font = Font(**font)
            else:
                font = Font()
        self['font'] = font
        self['select'] = select
        self['highlite'] = highlite
        self['dash'] = dash
        self['fixed'] = fixed
        self['edit'] = edit
        self['move'] = move
        self['rotate'] = rotate
        self['delete'] = delete
        self['include'] = include
        self['source'] = source
        self['text'] = text
        self['width'] = width

    @classmethod
    def none(cls):
        obj = cls()
        for k, v in obj.items():
            obj[k] = None
        return obj

    def format_value(self, k):
        if k in self:
            v = self[k]
            if k == 'background':
                if v:
                    return 'background'
            elif k == 'text':
                if v:
                    return f'text={{{v}}}'
            elif k == 'include':
                pass
            else:
                if isinstance(v, bool):
                    return f'{k}={int(v)}'
                elif v is None:
                    pass
                else:
                    return f'{k}={v}'

        return ''

    def __str__(self):
        out = []
        out = [self.format_value(k) for k in self.keys()]
        return ' '.join(filter(None, out))


if __name__ == '__main__':
    # x=Attributes()
    # y=str(x)
    #
    # print(y)
    z = Attributes.as_none()

    print(z)
