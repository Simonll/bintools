def joint_kwargs(**kwargs) -> str:
    return " ".join([k + " " + v for k, v in kwargs.items()])
