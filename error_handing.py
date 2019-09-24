#!/user/bin/env python
#coding=utf-8

from functools import wraps

class ErrorHandling(object):
    # @staticmethod
    pass 


class DecoratorHelper(object):
    @staticmethod
    def TRY_RETURN(except_retVal = False, exceptions = (Exception) ):
        def try_decotaator(fn):
            @wraps(fn)
            def decorated(*args, **kwargs):
                try:
                    return fn(*args, **kwargs)
                except exceptions as e:
                    print(e)
                    return except_retVal
            return  decorated
        return try_decotaator

    def TRY_RESUME_NEXT(func_or_var,default_value=None, exceptions = (Exception),is_ret_exception =False,*args, **kwargs):
            """
            :type default: object
            """
            if type(func_or_var) == type(lambda :x):
                func =  func_or_var
            else:
                func = lambda : func_or_var

            try:
                rtv = func(*args, **kwargs)
                if is_ret_exception:
                    return rtv,None
                else:
                    return rtv
            except exceptions as e:
                print(e)
                if is_ret_exception:
                    return default_value,e
                else:
                    return default_value


TRY_RESUME_NEXT = DecoratorHelper.TRY_RESUME_NEXT
TRY_RETURN = DecoratorHelper.TRY_RETURN

#
#  usage of DecoratorHelper
#
@DecoratorHelper.TRY_RETURN(except_retVal=None)
def fn(name= None):
    if name is None:
        raise Exception("name is none!")
    print("name is ",name)
    
    return {"name": name}


rtv = fn()
print("return: ",rtv)

rtv = fn("zhang1")
print("return: ",rtv)


#
#  usage of DecoratorHelper
#
# @DecoratorHelper.TRY_RETURN(except_retVal=None)
# def fn(name= None):
#     if name is None:
#         raise Exception("name is none!")
#     print("name is ",name)
#     return {"name": name}
#
# rtv = fn()
# print("return: ",rtv)
#
# rtv = fn("zhang1")
# print("return: ",rtv)