module_path = os.path.abspath(os.path.join('..', 'c0_preprocessing'))
if module_path not in sys.path:
    sys.path.append(module_path)
    
from c0 import connect_db