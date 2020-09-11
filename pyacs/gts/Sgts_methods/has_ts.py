###################################################################
def has_ts(self, code):
###################################################################
    """
    Tests whether a time series exists in the current Sgts instance
     
    :param code: 4-character code
    """
     
    def has_key(key, data):
        return True if key in data else False
    
    if has_key(code , self.__dict__):
        return True
    else:
        return False
