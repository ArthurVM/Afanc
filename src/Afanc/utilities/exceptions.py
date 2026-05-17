class AfancError(Exception):
    """Base class for exceptions in Afanc."""
    pass


class DirectoryExistsErrorAfanc(AfancError):
    """Raised when an output directory that should be new already exists."""
    def __init__(self, directory_path, message=None):
        self.directory_path = directory_path
        self.message = message or f"Output directory {directory_path} already exists! Please specify a new directory or delete the existing one."
        super().__init__(self.message)


class CommandExecutionError(AfancError):
    """Raised when an external command fails."""
    def __init__(self, command_str, return_code, stdout, stderr, message=None):
        self.command_str = command_str
        self.return_code = return_code
        self.stdout = stdout
        self.stderr = stderr
        self.message = message or f"Command '{command_str}' failed with return code {return_code}."
        super().__init__(self.message)

    def __str__(self):
        return f"{self.message}\nSTDOUT: {self.stdout}\nSTDERR: {self.stderr}"


class FileNotFoundErrorAfanc(AfancError, FileNotFoundError):
    """Raised when an essential file is not found."""
    def __init__(self, file_path, message=None):
        self.file_path = file_path
        self.message = message or f"Essential file not found: {file_path}"
        super().__init__(self.message)


class DirectoryNotFoundErrorAfanc(AfancError, FileNotFoundError):
    """Raised when an essential directory is not found."""
    def __init__(self, dir_path, message=None):
        self.dir_path = dir_path
        self.message = message or f"Essential directory not found: {dir_path}"
        super().__init__(self.message)


class InvalidFileFormatError(AfancError):
    """Raised when a file does not conform to the expected format."""
    def __init__(self, file_path, line_number=None, details=None, message=None):
        self.file_path = file_path
        self.line_number = line_number
        self.details = details
        _message = f"File '{file_path}' has an invalid format."
        if line_number is not None:
            _message += f" Error near line {line_number}."
        if details:
            _message += f" Details: {details}"
        self.message = message or _message
        super().__init__(self.message)

    
class SubtoolFailureError(AfancError):
    """Raised when a subtool fails to complete successfully."""
    def __init__(self, subtool_name, message=None):
        self.subtool_name = subtool_name
        self.message = message or f"Subtool '{subtool_name}' failed to complete."
        super().__init__(self.message)