# Product Guidelines

## Documentation & Tone
- **Concise & Practical:** Focus documentation and code comments on immediate application and "how-to." Technical explanations should be brief and aimed at helping the user achieve their goal quickly and correctly.

## Error Handling & Logging
- **Informative Logging:** Utilize robust logging (preferring `loguru`) to provide clear, actionable information about the script's execution. Errors must be clearly described, providing context to help the user resolve the issue.

## Code Organization & Design
- **Functional Decomposition:** Complex logic must be broken down into small, single-purpose functions. This ensures that the code is easier to understand, test, and maintain.

## Interface Design
- **Consistent CLI:** All scripts should provide a consistent and user-friendly command-line interface using a framework like `typer`. This includes helpful descriptions for arguments and options, ensuring a predictable experience for the user.
